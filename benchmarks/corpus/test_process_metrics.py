"""Independent stdlib tests: python -m unittest discover -s benchmarks/corpus -v."""

import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from typing import Any
import unittest
from unittest import mock

if __package__:
    from . import container_run, process_metrics
else:
    import container_run
    import process_metrics


class ProcessMetricsTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="process-metrics-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.stdout = self.root / "output" / "stdout.bin"
        self.stderr = self.root / "errors" / "stderr.bin"

    def run_child(self, source, timeout=5.0):
        return process_metrics.run_process(
            [sys.executable, "-c", source], self.root,
            self.stdout, self.stderr, timeout,
        )

    def assert_metrics(self, result, status, returncode):
        self.assertEqual(result["status"], status, result)
        self.assertEqual(result["returncode"], returncode, result)
        self.assertIs(type(result["returncode"]), int)
        self.assertIs(type(result["wall_seconds"]), float)
        self.assertTrue(math.isfinite(result["wall_seconds"]))
        self.assertGreaterEqual(result["wall_seconds"], 0.0)
        self.assertLess(result["wall_seconds"], 10.0)
        self.assertIsNone(result["cleanup_error"], result)
        self.assertIsInstance(result["cleanup_method"], str)
        cpu = result["cpu_seconds"]
        if cpu is None:
            self.assertIsInstance(result["cpu_error"], str)
            self.assertTrue(result["cpu_error"])
        else:
            self.assertIs(type(cpu), float)
            self.assertTrue(math.isfinite(cpu))
            self.assertGreaterEqual(cpu, 0.0)
            # All test children are single-threaded: seconds, not 100ns ticks.
            self.assertLessEqual(cpu, result["wall_seconds"] + 0.25)
            self.assertIsInstance(result["cpu_method"], str)
            self.assertIsNone(result["cpu_error"])
        rss = result["peak_rss_bytes"]
        if rss is None:
            self.assertIsInstance(result["memory_error"], str)
            self.assertTrue(result["memory_error"])
        else:
            self.assertIs(type(rss), int)
            self.assertGreater(rss, 1024 * 1024)
            self.assertIsInstance(result["memory_method"], str)
            self.assertIsNone(result["memory_error"])
        if os.name == "nt":
            self.assertIsNotNone(cpu, result)
            self.assertIsNotNone(rss, result)
            self.assertIn("process handle", result["cpu_method"])
            self.assertIn("PeakWorkingSetSize", result["memory_method"])

    def assert_process_stopped(self, pid):
        """Check without ever terminating a PID that might have been reused."""
        if os.name == "nt":
            import ctypes
            from ctypes import wintypes

            kernel = ctypes.WinDLL("kernel32", use_last_error=True)
            kernel.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
            kernel.OpenProcess.restype = wintypes.HANDLE
            kernel.WaitForSingleObject.argtypes = [wintypes.HANDLE, wintypes.DWORD]
            kernel.WaitForSingleObject.restype = wintypes.DWORD
            kernel.CloseHandle.argtypes = [wintypes.HANDLE]
            kernel.CloseHandle.restype = wintypes.BOOL
            handle = kernel.OpenProcess(0x00100000, False, pid)  # SYNCHRONIZE
            if not handle:
                self.assertEqual(ctypes.get_last_error(), 87)  # ERROR_INVALID_PARAMETER: gone
                return
            try:
                self.assertEqual(kernel.WaitForSingleObject(handle, 1000), 0, f"PID {pid} alive")
            finally:
                kernel.CloseHandle(handle)
        elif sys.platform.startswith("linux"):
            deadline = time.perf_counter() + 1.0
            while True:
                try:
                    # Orphaned zombies are exited, though their init may not reap promptly.
                    stat = Path(f"/proc/{pid}/stat").read_text()
                except FileNotFoundError:
                    return
                if stat.rsplit(")", 1)[1].split()[0] == "Z":
                    return
                if time.perf_counter() >= deadline:
                    self.fail(f"PID {pid} is still running")
                time.sleep(0.01)
        else:
            self.skipTest("Reliable non-child exit observation implemented for Windows/Linux only")

    def test_successful_empty_and_very_short_children(self):
        for _ in range(5):
            with self.subTest():
                result = self.run_child("pass")
                self.assert_metrics(result, "success", 0)
                self.assertIsNone(result["error"])
                self.assertEqual(self.stdout.read_bytes(), b"")
                self.assertEqual(self.stderr.read_bytes(), b"")

    def test_nonzero_exit_and_binary_output(self):
        result = self.run_child(
            "import sys; sys.stdout.buffer.write(b'out\\x00\\xff\\n'); "
            "sys.stderr.buffer.write(b'err\\x00\\xfe\\n'); sys.exit(7)"
        )
        self.assert_metrics(result, "error", 7)
        self.assertEqual(self.stdout.read_bytes(), b"out\x00\xff\n")
        self.assertEqual(self.stderr.read_bytes(), b"err\x00\xfe\n")

    def test_large_output_goes_to_files_without_pipe_deadlock(self):
        result = self.run_child(
            "import sys; sys.stdout.buffer.write(b'o' * (2 * 1024 * 1024)); "
            "sys.stderr.buffer.write(b'e' * (2 * 1024 * 1024))"
        )
        self.assert_metrics(result, "success", 0)
        self.assertEqual(self.stdout.stat().st_size, 2 * 1024 * 1024)
        self.assertEqual(self.stderr.stat().st_size, 2 * 1024 * 1024)

    def test_memory_allocation_is_in_bytes_and_per_child(self):
        allocated = self.run_child(
            "data = bytearray(64 * 1024 * 1024); "
            "data[::4096] = b'x' * (len(data) // 4096)"
        )
        self.assert_metrics(allocated, "success", 0)
        small = self.run_child("pass")
        self.assert_metrics(small, "success", 0)
        if allocated["peak_rss_bytes"] is not None:
            self.assertGreaterEqual(allocated["peak_rss_bytes"], 64 * 1024 * 1024)
            self.assertLess(small["peak_rss_bytes"], allocated["peak_rss_bytes"] - 32 * 1024 * 1024)

    def test_cpu_seconds_are_for_child_not_measuring_process(self):
        # Parent CPU must not be substituted for a short child's usage.
        parent_deadline = time.perf_counter() + 0.15
        while time.perf_counter() < parent_deadline:
            pass
        busy = self.run_child(
            "import time; end = time.perf_counter() + 0.2\n"
            "while time.perf_counter() < end: pass"
        )
        self.assert_metrics(busy, "success", 0)
        idle = self.run_child("pass")
        self.assert_metrics(idle, "success", 0)
        if busy["cpu_seconds"] is not None:
            self.assertGreater(busy["cpu_seconds"], 0.02)
            self.assertLess(idle["cpu_seconds"], busy["cpu_seconds"])

    def test_timeout_captures_output_and_reaps_child(self):
        result = self.run_child(
            "import os,sys,time; print(os.getpid(), flush=True); "
            "print('before timeout', file=sys.stderr, flush=True); time.sleep(10)",
            timeout=1.5,
        )
        self.assertEqual(result["status"], "timeout", result)
        self.assertIs(type(result["returncode"]), int)
        self.assertNotEqual(result["returncode"], 0)
        self.assert_metrics(result, "timeout", result["returncode"])
        self.assertGreaterEqual(result["wall_seconds"], 1.4)
        self.assertLess(result["wall_seconds"], 4.0)
        self.assertIn("timeout", result["error"])
        self.assertEqual(self.stderr.read_bytes(), b"before timeout\r\n" if os.name == "nt" else b"before timeout\n")
        self.assert_process_stopped(int(self.stdout.read_text().strip()))

    @unittest.skipUnless(os.name == "nt" or sys.platform.startswith("linux"), "Exit check needs Windows/Linux")
    def test_timeout_terminates_descendant_without_harming_unrelated_child(self):
        # A control process is deliberately outside the runner's job/session.
        control = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(10)"],
            stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        try:
            result = self.run_child(
                "import os,subprocess,sys,time; "
                "child = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(10)']); "
                "print(os.getpid(), child.pid, flush=True); time.sleep(10)",
                timeout=1.5,
            )
            self.assert_metrics(result, "timeout", result["returncode"])
            pids = [int(value) for value in self.stdout.read_text().split()]
            self.assertEqual(len(pids), 2)
            for pid in pids:
                self.assert_process_stopped(pid)
            self.assertIsNone(control.poll(), "Unrelated process was terminated")
        finally:
            if control.poll() is None:
                control.kill()
            control.wait(timeout=3)

    @unittest.skipUnless(os.name == "nt", "Job Objects are Windows-specific")
    def test_windows_normal_exit_cleans_up_descendant(self):
        result = self.run_child(
            "import subprocess,sys; "
            "child = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(10)']); "
            "print(child.pid, flush=True)"
        )
        self.assert_metrics(result, "success", 0)
        self.assert_process_stopped(int(self.stdout.read_text().strip()))

    def test_startup_is_part_of_timeout(self):
        result = self.run_child("pass", timeout=1e-9)
        self.assert_metrics(result, "timeout", result["returncode"])
        self.assertGreater(result["wall_seconds"], 1e-9)

    def test_argv_quoting_cwd_and_stdin_eof(self):
        arguments = ["", "two words", 'quote"inside', "backslash\\", "snowman \u2603"]
        result = process_metrics.run_process(
            [sys.executable, "-c",
             "import os,sys; sys.stdout.reconfigure(encoding='utf-8'); "
             "print(repr(sys.argv[1:])); print(os.getcwd()); "
             "print(len(sys.stdin.buffer.read()))", *arguments],
            self.root, self.stdout, self.stderr, 5.0,
        )
        self.assert_metrics(result, "success", 0)
        lines = self.stdout.read_text(encoding="utf-8").splitlines()
        self.assertEqual(lines, [repr(arguments), str(self.root), "0"])

    def test_large_finite_timeout_does_not_overflow_platform_wait(self):
        result = self.run_child("pass", timeout=1e308)
        self.assert_metrics(result, "success", 0)

    def test_invalid_timeouts_have_no_output_side_effects(self):
        for value in (0, -1, float("inf"), float("-inf"), float("nan"), True):
            with self.subTest(timeout=value), self.assertRaises(ValueError):
                self.run_child("pass", timeout=value)
        self.assertEqual(list(self.root.iterdir()), [])

    def test_stdout_and_stderr_cannot_be_the_same_path(self):
        for path in (self.stdout, self.stdout.parent / "sub" / ".." / self.stdout.name):
            with self.subTest(path=path), self.assertRaises(ValueError):
                process_metrics.run_process(
                    [sys.executable, "-c", "pass"], self.root,
                    self.stdout, path, 5.0,
                )
        self.assertEqual(list(self.root.iterdir()), [])

    def test_existing_hardlinked_outputs_are_rejected_without_truncation(self):
        first, second = self.root / "first", self.root / "second"
        first.write_bytes(b"preserve me")
        try:
            os.link(first, second)
        except OSError as exc:
            self.skipTest(f"Hard links unavailable: {exc}")
        with self.assertRaises(ValueError):
            process_metrics.run_process(
                [sys.executable, "-c", "pass"], self.root, first, second, 5.0,
            )
        self.assertEqual(first.read_bytes(), b"preserve me")

    def test_invalid_argv_has_no_output_side_effects(self):
        invalid_cases: tuple[Any, ...] = (
            [], [""], [sys.executable, "bad\0argument"], [123], "not a list",
        )
        for argv in invalid_cases:
            with self.subTest(argv=argv), self.assertRaises(ValueError):
                process_metrics.run_process(argv, self.root, self.stdout, self.stderr, 5.0)
        self.assertEqual(list(self.root.iterdir()), [])

    def test_launch_error_reports_unknown_metrics_and_only_makes_output_parents(self):
        result = process_metrics.run_process(
            [str(self.root / "missing-executable")], self.root,
            self.stdout, self.stderr, 5.0,
        )
        self.assertEqual(result["status"], "error")
        self.assertIsNone(result["returncode"])
        self.assertIsNone(result["cpu_seconds"])
        self.assertIsNone(result["peak_rss_bytes"])
        self.assertTrue(result["error"])
        self.assertTrue(result["cpu_error"])
        self.assertTrue(result["memory_error"])
        self.assertEqual({p.name for p in self.root.iterdir()}, {"output", "errors"})
        self.assertEqual(self.stdout.read_bytes(), b"")
        self.assertEqual(self.stderr.read_bytes(), b"")

    def test_output_open_error_reports_measured_wall_and_unknown_child_metrics(self):
        self.stderr.mkdir(parents=True)
        result = self.run_child("pass")
        self.assertEqual(result["status"], "error")
        self.assertIsNone(result["returncode"])
        self.assertTrue(result["error"])
        self.assertIs(type(result["wall_seconds"]), float)
        self.assertGreater(result["wall_seconds"], 0.0)
        self.assertIsNone(result["cpu_seconds"])
        self.assertIsNone(result["peak_rss_bytes"])

    def test_missing_cwd_is_not_created(self):
        missing = self.root / "not-a-working-directory"
        result = process_metrics.run_process(
            [sys.executable, "-c", "pass"], missing, self.stdout, self.stderr, 5.0,
        )
        self.assertEqual(result["status"], "error")
        self.assertFalse(missing.exists())
        self.assertIsNone(result["returncode"])

    def test_existing_outputs_are_truncated(self):
        self.run_child("print('old output')")
        result = self.run_child("pass")
        self.assert_metrics(result, "success", 0)
        self.assertEqual(self.stdout.read_bytes(), b"")
        self.assertEqual(self.stderr.read_bytes(), b"")

    @unittest.skipUnless(os.name == "nt", "Windows API failure injection")
    def test_assignment_failure_never_resumes_child(self):
        import ctypes

        marker = self.root / "child-ran"
        api = process_metrics._WindowsAPI()

        def deny_assignment(*args):
            ctypes.set_last_error(5)
            return False

        setattr(api, "AssignProcessToJobObject", deny_assignment)
        with mock.patch.object(process_metrics, "_WindowsAPI", return_value=api):
            result = self.run_child(f"from pathlib import Path; Path({str(marker)!r}).touch()")
        self.assertEqual(result["status"], "error")
        self.assertIn("AssignProcessToJobObject", result["error"])
        self.assertIn("not resumed", result["error"])
        self.assertIsNone(result["cleanup_error"])
        self.assertFalse(marker.exists())

    @unittest.skipUnless(os.name == "nt", "Windows API failure injection")
    def test_metric_query_failures_report_null_not_zero(self):
        import ctypes

        for function, metric, error_key in (
            ("GetProcessMemoryInfo", "peak_rss_bytes", "memory_error"),
            ("GetProcessTimes", "cpu_seconds", "cpu_error"),
        ):
            with self.subTest(function=function):
                api = process_metrics._WindowsAPI()

                def fail_query(*args):
                    ctypes.set_last_error(5)
                    return False

                setattr(api, function, fail_query)
                with mock.patch.object(process_metrics, "_WindowsAPI", return_value=api):
                    result = self.run_child("pass")
                self.assertEqual(result["status"], "success", result)
                self.assertEqual(result["returncode"], 0)
                self.assertIsNone(result[metric])
                self.assertIn(function, result[error_key])


class EvidenceTests(unittest.TestCase):
    """Failure classification and cgroup evidence never fabricate a zero."""

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="evidence-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def test_tail_text_reports_full_size_and_bounded_tail(self):
        self.assertEqual(process_metrics.tail_text(self.root / "missing"), ("", 0))
        small = self.root / "small"
        small.write_bytes(b"hello")
        self.assertEqual(process_metrics.tail_text(small), ("hello", 5))
        large = self.root / "large"
        large.write_bytes(b"a" * 100 + b"tail")
        self.assertEqual(process_metrics.tail_text(large, max_bytes=4), ("tail", 104))
        binary = self.root / "binary"
        binary.write_bytes(b"ok\xff")
        self.assertEqual(process_metrics.tail_text(binary), ("ok\ufffd", 3))

    def test_classify_exit_requires_evidence_for_oom(self):
        cases = (
            (0, "", None, "success"),
            (1, "", None, "nonzero"),
            (1, "java.lang.OutOfMemoryError: Java heap space", None, "oom-heap"),
            (137, "", None, "killed"),
            (137, "", 0, "killed"),
            (137, "", 1, "oom-killed"),
            (1, "", 1, "oom-killed"),
            (-9, "", None, "killed"),
            (9, "", None, "killed"),
            (1, "", 0, "nonzero"),
            (None, "", None, "nonzero"),
        )
        for returncode, stderr, kills, expected in cases:
            with self.subTest(returncode=returncode, stderr=stderr, kills=kills):
                self.assertEqual(
                    process_metrics.classify_exit(returncode, stderr, kills),
                    expected,
                )

    def test_failure_messages_name_the_evidence(self):
        self.assertIn(
            "OOM kill", container_run.failure_message("oom-killed", -9)
        )
        self.assertIn(
            "OutOfMemoryError", container_run.failure_message("oom-heap", 1)
        )
        self.assertIn(
            "without cgroup OOM evidence",
            container_run.failure_message("killed", 137),
        )
        self.assertIn(
            "status 7", container_run.failure_message("nonzero", 7)
        )

    def test_read_cgroup_memory_v2(self):
        (self.root / "memory.peak").write_text("123456\n", encoding="ascii")
        (self.root / "memory.events").write_text(
            "low 0\nhigh 1\nmax 2\noom 3\noom_kill 4\n", encoding="ascii"
        )
        self.assertEqual(process_metrics.read_cgroup_memory(self.root), {
            "cgroup_memory_peak_bytes": 123456,
            "cgroup_oom_kills": 4,
            "cgroup_memory_method": "cgroup v2 memory.peak/memory.events",
        })

    def test_read_cgroup_memory_v2_without_events_keeps_peak(self):
        (self.root / "memory.peak").write_text("99", encoding="ascii")
        result = process_metrics.read_cgroup_memory(self.root)
        self.assertEqual(result["cgroup_memory_peak_bytes"], 99)
        self.assertIsNone(result["cgroup_oom_kills"])
        (self.root / "memory.peak").write_text("not a number", encoding="ascii")
        result = process_metrics.read_cgroup_memory(self.root)
        self.assertIsNone(result["cgroup_memory_peak_bytes"])

    def test_read_cgroup_memory_v1_reports_peak_but_not_failcnt_as_kills(self):
        memory = self.root / "memory"
        memory.mkdir()
        (memory / "memory.max_usage_in_bytes").write_text("777\n", encoding="ascii")
        (memory / "memory.failcnt").write_text("12\n", encoding="ascii")
        result = process_metrics.read_cgroup_memory(self.root)
        self.assertEqual(result["cgroup_memory_peak_bytes"], 777)
        self.assertIsNone(result["cgroup_oom_kills"])
        self.assertIn("v1", result["cgroup_memory_method"])

    def test_read_cgroup_memory_absent_is_unknown(self):
        empty = self.root / "empty"
        empty.mkdir()
        result = process_metrics.read_cgroup_memory(empty)
        self.assertEqual(result, {
            "cgroup_memory_peak_bytes": None,
            "cgroup_oom_kills": None,
            "cgroup_memory_method": None,
        })


if __name__ == "__main__":
    unittest.main()
