"""Disk-only output capture and per-process CLI measurements (standard library).

Windows children start suspended and enter a kill-on-close Job Object before
running. Descendants cannot break away from that job. Job setup/assignment
failure is a launch error, not permission to run an uncontained child. Metrics
come from the original process handle, queried after exit and before closing it;
CPU and peak working set describe the direct child, not the entire job.

POSIX uses a new session and kills its initial process group on timeout, before
reaping the leader (avoiding PID-reuse races). Descendants that change groups or
sessions are not covered, and descendants are not tracked after a normal exit.
Where available, wait4 supplies per-child rusage, which can include descendants
that child reaped. No cumulative RUSAGE_CHILDREN measurements are used.

Invalid arguments raise ValueError. Launch/measurement/cleanup failures are
reported in the result. Unknown metrics are None, never synthetic zeroes.
Output paths are relative to the caller, not cwd. Output files are truncated.
The wall interval includes process startup, waiting, and bounded cleanup, but
not argument validation or opening output files. If output setup fails, wall
instead records that failed setup attempt. Cleanup has a separate two-second
budget; a timeout is therefore not a strict bound on total wall time.
"""

from __future__ import annotations

import ctypes
from ctypes import wintypes
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


_CLEANUP_SECONDS = 2.0


def run_process(
    argv: list[str],
    cwd: Path,
    stdout_path: Path,
    stderr_path: Path,
    timeout_seconds: float,
) -> dict:
    """Run one CLI, returning status, timings, peak RSS, and measurement errors.

    returncode is an int if the child was reaped/observed exited, otherwise None.
    status is success only for exit code zero and successful cleanup. Unavailable
    metrics do not change the exit status; see memory_error and cpu_error.
    CPU seconds and RSS bytes are floats and integers respectively when known.
    """
    if (
        isinstance(timeout_seconds, bool)
        or not math.isfinite(timeout_seconds)
        or timeout_seconds <= 0
    ):
        raise ValueError("timeout_seconds must be finite and positive")
    if (
        not isinstance(argv, list)
        or not argv
        or not argv[0]
        or any(not isinstance(arg, str) or "\0" in arg for arg in argv)
    ):
        raise ValueError("argv must be a nonempty list of NUL-free strings")
    cwd, stdout_path, stderr_path = Path(cwd), Path(stdout_path), Path(stderr_path)
    if "\0" in str(cwd):
        raise ValueError("cwd must not contain NUL")
    if stdout_path.resolve() == stderr_path.resolve():
        raise ValueError("stdout_path and stderr_path must be different files")
    try:
        same_file = stdout_path.samefile(stderr_path)
    except FileNotFoundError:
        same_file = False
    if same_file:
        raise ValueError("stdout_path and stderr_path must be different files")

    result = {
        "returncode": None,
        "status": "error",
        "wall_seconds": 0.0,
        "cpu_seconds": None,
        "cpu_method": None,
        "cpu_error": "Child CPU usage was not measured (process did not complete).",
        "peak_rss_bytes": None,
        "memory_method": None,
        "memory_error": "Child peak RSS was not measured (process did not complete).",
        "error": None,
        "cleanup_method": None,
        "cleanup_error": None,
    }
    started = time.perf_counter()
    try:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
        stderr_path.parent.mkdir(parents=True, exist_ok=True)
        with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
            started = time.perf_counter()
            deadline = started + timeout_seconds
            if os.name == "nt":
                _run_windows(argv, cwd, stdout, stderr, deadline, result)
            elif os.name == "posix":
                _run_posix(argv, cwd, stdout, stderr, deadline, result)
            else:
                result["error"] = f"Unsupported process platform: {os.name}"
    except OSError as exc:
        result["status"] = "error"
        result["error"] = str(exc)
    finally:
        result["wall_seconds"] = float(time.perf_counter() - started)
    return result


def _record_cleanup_error(result, message):
    previous = result["cleanup_error"]
    result["cleanup_error"] = f"{previous}; {message}" if previous else message


def _run_posix(argv, cwd, stdout, stderr, deadline, result):
    result["cleanup_method"] = (
        "On timeout: SIGKILL to the owned initial process group before reaping; "
        "descendants leaving that group and normal-exit descendants are not covered."
    )
    with open(os.devnull, "rb") as stdin:
        process = subprocess.Popen(
            argv, cwd=cwd, stdin=stdin, stdout=stdout, stderr=stderr,
            start_new_session=True, close_fds=True,
        )
    usage = None
    reaped = False

    def reap():
        nonlocal usage, reaped
        if hasattr(os, "wait4"):
            pid, status, child_usage = os.wait4(process.pid, os.WNOHANG)
            if pid:
                process.returncode = os.waitstatus_to_exitcode(status)
                usage = child_usage
                reaped = True
        else:
            reaped = process.poll() is not None
        return reaped

    try:
        if time.perf_counter() >= deadline:
            result["status"] = "timeout"
            result["error"] = "Process startup exceeded timeout_seconds."
        while result["status"] != "timeout" and not reap():
            remaining = deadline - time.perf_counter()
            if remaining <= 0:
                result["status"] = "timeout"
                result["error"] = "Process exceeded timeout_seconds."
                break
            time.sleep(min(0.01, remaining))
    finally:
        if not reaped:
            # The unreaped session leader reserves its PID/process-group ID.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            except OSError as exc:
                _record_cleanup_error(result, f"killpg: {exc}")
                try:
                    process.kill()
                except ProcessLookupError:
                    pass
                except OSError as kill_exc:
                    _record_cleanup_error(result, f"kill child: {kill_exc}")
            cleanup_deadline = time.perf_counter() + _CLEANUP_SECONDS
            while not reap() and time.perf_counter() < cleanup_deadline:
                time.sleep(0.01)
            if not reaped:
                _record_cleanup_error(result, "Child did not exit within cleanup budget.")

    result["returncode"] = process.returncode
    if result["status"] != "timeout":
        result["status"] = (
            "success" if process.returncode == 0 and not result["cleanup_error"] else "error"
        )
    if usage is None:
        reason = "Per-child os.wait4 rusage is unavailable or the child was not reaped."
        result["cpu_error"] = reason
        result["memory_error"] = reason
        return
    result["cpu_seconds"] = float(usage.ru_utime + usage.ru_stime)
    result["cpu_method"] = "os.wait4 rusage: child, including descendants it reaped"
    result["cpu_error"] = None
    if sys.platform == "darwin":
        multiplier = 1
    elif sys.platform.startswith("linux"):
        multiplier = 1024
    else:
        result["memory_error"] = f"Unknown wait4 ru_maxrss units on {sys.platform}."
        return
    peak = int(usage.ru_maxrss) * multiplier
    if peak > 0:
        result["peak_rss_bytes"] = peak
        result["memory_method"] = (
            "os.wait4 ru_maxrss: child, including descendants it reaped; converted to bytes"
        )
        result["memory_error"] = None
    else:
        result["memory_error"] = "os.wait4 reported no usable peak RSS."


if sys.platform == "win32":
    class _IO_COUNTERS(ctypes.Structure):
        _fields_ = [(name, ctypes.c_ulonglong) for name in (
            "ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
            "ReadTransferCount", "WriteTransferCount", "OtherTransferCount",
        )]

    class _BASIC_LIMITS(ctypes.Structure):
        _fields_ = [
            ("PerProcessUserTimeLimit", ctypes.c_longlong),
            ("PerJobUserTimeLimit", ctypes.c_longlong),
            ("LimitFlags", wintypes.DWORD),
            ("MinimumWorkingSetSize", ctypes.c_size_t),
            ("MaximumWorkingSetSize", ctypes.c_size_t),
            ("ActiveProcessLimit", wintypes.DWORD),
            ("Affinity", ctypes.c_size_t),
            ("PriorityClass", wintypes.DWORD),
            ("SchedulingClass", wintypes.DWORD),
        ]

    class _EXTENDED_LIMITS(ctypes.Structure):
        _fields_ = [
            ("BasicLimitInformation", _BASIC_LIMITS), ("IoInfo", _IO_COUNTERS),
            ("ProcessMemoryLimit", ctypes.c_size_t), ("JobMemoryLimit", ctypes.c_size_t),
            ("PeakProcessMemoryUsed", ctypes.c_size_t), ("PeakJobMemoryUsed", ctypes.c_size_t),
        ]

    class _JOB_ACCOUNTING(ctypes.Structure):
        _fields_ = [
            ("TotalUserTime", ctypes.c_longlong), ("TotalKernelTime", ctypes.c_longlong),
            ("ThisPeriodTotalUserTime", ctypes.c_longlong),
            ("ThisPeriodTotalKernelTime", ctypes.c_longlong),
            ("TotalPageFaultCount", wintypes.DWORD), ("TotalProcesses", wintypes.DWORD),
            ("ActiveProcesses", wintypes.DWORD), ("TotalTerminatedProcesses", wintypes.DWORD),
        ]

    class _MEMORY_COUNTERS(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD)] + [
            (name, ctypes.c_size_t) for name in (
                "PeakWorkingSetSize", "WorkingSetSize", "QuotaPeakPagedPoolUsage",
                "QuotaPagedPoolUsage", "QuotaPeakNonPagedPoolUsage",
                "QuotaNonPagedPoolUsage", "PagefileUsage", "PeakPagefileUsage",
            )
        ]

    class _STARTUPINFO(ctypes.Structure):
        _fields_ = [
            ("cb", wintypes.DWORD), ("lpReserved", wintypes.LPWSTR),
            ("lpDesktop", wintypes.LPWSTR), ("lpTitle", wintypes.LPWSTR),
            ("dwX", wintypes.DWORD), ("dwY", wintypes.DWORD),
            ("dwXSize", wintypes.DWORD), ("dwYSize", wintypes.DWORD),
            ("dwXCountChars", wintypes.DWORD), ("dwYCountChars", wintypes.DWORD),
            ("dwFillAttribute", wintypes.DWORD), ("dwFlags", wintypes.DWORD),
            ("wShowWindow", wintypes.WORD), ("cbReserved2", wintypes.WORD),
            ("lpReserved2", ctypes.POINTER(ctypes.c_ubyte)),
            ("hStdInput", wintypes.HANDLE), ("hStdOutput", wintypes.HANDLE),
            ("hStdError", wintypes.HANDLE),
        ]

    class _STARTUPINFOEX(ctypes.Structure):
        _fields_ = [("StartupInfo", _STARTUPINFO), ("lpAttributeList", ctypes.c_void_p)]

    class _PROCESS_INFORMATION(ctypes.Structure):
        _fields_ = [
            ("hProcess", wintypes.HANDLE), ("hThread", wintypes.HANDLE),
            ("dwProcessId", wintypes.DWORD), ("dwThreadId", wintypes.DWORD),
        ]

    class _WindowsAPI:
        def __init__(self):
            kernel = ctypes.WinDLL("kernel32", use_last_error=True)
            psapi = ctypes.WinDLL("psapi", use_last_error=True)
            H, D, P, B = wintypes.HANDLE, wintypes.DWORD, ctypes.c_void_p, wintypes.BOOL

            def bind(name, restype, *argtypes):
                function = getattr(kernel, name)
                function.restype, function.argtypes = restype, argtypes
                return function

            self.CreateJobObjectW = bind("CreateJobObjectW", H, P, wintypes.LPCWSTR)
            self.SetInformationJobObject = bind("SetInformationJobObject", B, H, ctypes.c_int, P, D)
            self.AssignProcessToJobObject = bind("AssignProcessToJobObject", B, H, H)
            self.TerminateJobObject = bind("TerminateJobObject", B, H, wintypes.UINT)
            self.QueryInformationJobObject = bind("QueryInformationJobObject", B, H, ctypes.c_int, P, D, P)
            self.InitializeProcThreadAttributeList = bind("InitializeProcThreadAttributeList", B, P, D, D, P)
            self.UpdateProcThreadAttribute = bind(
                "UpdateProcThreadAttribute", B, P, D, ctypes.c_size_t, P, ctypes.c_size_t, P, P,
            )
            self.DeleteProcThreadAttributeList = bind("DeleteProcThreadAttributeList", None, P)
            self.GetCurrentProcess = bind("GetCurrentProcess", H)
            self.DuplicateHandle = bind("DuplicateHandle", B, H, H, H, P, D, B, D)
            self.CreateProcessW = bind(
                "CreateProcessW", B, wintypes.LPCWSTR, wintypes.LPWSTR, P, P, B, D, P,
                wintypes.LPCWSTR, P, P,
            )
            self.ResumeThread = bind("ResumeThread", D, H)
            self.TerminateProcess = bind("TerminateProcess", B, H, wintypes.UINT)
            self.WaitForSingleObject = bind("WaitForSingleObject", D, H, D)
            self.GetExitCodeProcess = bind("GetExitCodeProcess", B, H, P)
            self.GetProcessTimes = bind("GetProcessTimes", B, H, P, P, P, P)
            self.CloseHandle = bind("CloseHandle", B, H)
            self.GetProcessMemoryInfo = psapi.GetProcessMemoryInfo
            self.GetProcessMemoryInfo.restype = B
            self.GetProcessMemoryInfo.argtypes = [H, P, D]


def _windows_error(operation):
    return f"{operation}: {ctypes.WinError(ctypes.get_last_error())}"


def _create_suspended(api, argv, cwd, stdout, stderr):
    """Use an explicit inheritance list; keep the primary thread handle to resume."""
    import msvcrt

    size = ctypes.c_size_t()
    api.InitializeProcThreadAttributeList(None, 1, 0, ctypes.byref(size))
    if not size.value:
        raise OSError(_windows_error("InitializeProcThreadAttributeList(size)"))
    attributes = ctypes.create_string_buffer(size.value)
    if not api.InitializeProcThreadAttributeList(attributes, 1, 0, ctypes.byref(size)):
        raise OSError(_windows_error("InitializeProcThreadAttributeList"))
    inherited = []
    try:
        with open(os.devnull, "rb") as stdin:
            current = api.GetCurrentProcess()
            for stream in (stdin, stdout, stderr):
                duplicate = wintypes.HANDLE()
                if not api.DuplicateHandle(
                    current, msvcrt.get_osfhandle(stream.fileno()), current,
                    ctypes.byref(duplicate), 0, True, 2,  # DUPLICATE_SAME_ACCESS
                ):
                    raise OSError(_windows_error("DuplicateHandle"))
                inherited.append(duplicate.value)
            handles = (wintypes.HANDLE * len(inherited))(*inherited)
            if not api.UpdateProcThreadAttribute(
                attributes, 0, 0x00020002, handles, ctypes.sizeof(handles), None, None,
            ):  # PROC_THREAD_ATTRIBUTE_HANDLE_LIST
                raise OSError(_windows_error("UpdateProcThreadAttribute(handle list)"))
            startup = _STARTUPINFOEX()
            startup.StartupInfo.cb = ctypes.sizeof(startup)
            startup.StartupInfo.dwFlags = 0x100  # STARTF_USESTDHANDLES
            (startup.StartupInfo.hStdInput, startup.StartupInfo.hStdOutput,
             startup.StartupInfo.hStdError) = inherited
            startup.lpAttributeList = ctypes.cast(attributes, ctypes.c_void_p)
            info = _PROCESS_INFORMATION()
            command = ctypes.create_unicode_buffer(subprocess.list2cmdline(argv))
            # CREATE_SUSPENDED | EXTENDED_STARTUPINFO_PRESENT | CREATE_NO_WINDOW
            if not api.CreateProcessW(
                None, command, None, None, True, 0x00000004 | 0x00080000 | 0x08000000,
                None, str(cwd), ctypes.byref(startup), ctypes.byref(info),
            ):
                raise OSError(_windows_error("CreateProcessW"))
            return info
    finally:
        for handle in inherited:
            api.CloseHandle(handle)
        api.DeleteProcThreadAttributeList(attributes)


def _wait_windows(api, handle, deadline):
    while True:
        remaining = max(0.0, deadline - time.perf_counter())
        milliseconds = math.ceil(min(remaining, 0xFFFFFFFE / 1000) * 1000)
        state = api.WaitForSingleObject(handle, milliseconds)
        if state == 0:  # WAIT_OBJECT_0
            return True
        if state != 258:  # WAIT_TIMEOUT
            raise OSError(_windows_error("WaitForSingleObject"))
        if time.perf_counter() >= deadline:
            return False


def _cleanup_windows(api, job, info, assigned, result):
    deadline = time.perf_counter() + _CLEANUP_SECONDS
    if assigned:
        if not api.TerminateJobObject(job, 1):
            _record_cleanup_error(result, _windows_error("TerminateJobObject"))
    elif info is not None:
        # Assignment/resume never succeeded: this child is still suspended.
        if not api.TerminateProcess(info.hProcess, 1):
            _record_cleanup_error(result, _windows_error("TerminateProcess(suspended child)"))
    exited = False
    if info is not None:
        try:
            exited = _wait_windows(api, info.hProcess, deadline)
            if not exited:
                _record_cleanup_error(result, "Child did not exit within cleanup budget.")
        except OSError as exc:
            _record_cleanup_error(result, str(exc))
    if assigned:
        accounting = _JOB_ACCOUNTING()
        while True:
            if not api.QueryInformationJobObject(
                job, 1, ctypes.byref(accounting), ctypes.sizeof(accounting), None,
            ):  # JobObjectBasicAccountingInformation
                _record_cleanup_error(result, _windows_error("QueryInformationJobObject"))
                break
            if accounting.ActiveProcesses == 0:
                break
            if time.perf_counter() >= deadline:
                _record_cleanup_error(result, "Job still has active processes after cleanup budget.")
                break
            time.sleep(min(0.01, max(0.0, deadline - time.perf_counter())))
    return exited


def _measure_windows(api, process, result):
    code = wintypes.DWORD()
    if api.GetExitCodeProcess(process, ctypes.byref(code)):
        result["returncode"] = int(code.value)
    else:
        result["error"] = _windows_error("GetExitCodeProcess")
    creation, exit_time, kernel, user = (wintypes.FILETIME() for _ in range(4))
    if api.GetProcessTimes(
        process, ctypes.byref(creation), ctypes.byref(exit_time),
        ctypes.byref(kernel), ctypes.byref(user),
    ):
        ticks = sum((value.dwHighDateTime << 32) | value.dwLowDateTime for value in (kernel, user))
        result["cpu_seconds"] = float(ticks / 10_000_000)
        result["cpu_method"] = "GetProcessTimes kernel + user, direct child process handle"
        result["cpu_error"] = None
    else:
        result["cpu_error"] = _windows_error("GetProcessTimes")
    counters = _MEMORY_COUNTERS()
    counters.cb = ctypes.sizeof(counters)
    if api.GetProcessMemoryInfo(process, ctypes.byref(counters), ctypes.sizeof(counters)):
        if counters.PeakWorkingSetSize > 0:
            result["peak_rss_bytes"] = int(counters.PeakWorkingSetSize)
            result["memory_method"] = (
                "GetProcessMemoryInfo PROCESS_MEMORY_COUNTERS.PeakWorkingSetSize, "
                "direct child process handle; bytes"
            )
            result["memory_error"] = None
        else:
            result["memory_error"] = "GetProcessMemoryInfo reported no usable peak working set."
    else:
        result["memory_error"] = _windows_error("GetProcessMemoryInfo")


def _run_windows(argv, cwd, stdout, stderr, deadline, result):
    result["cleanup_method"] = (
        "Suspended child assigned before resume to a KILL_ON_JOB_CLOSE Job Object; "
        "TerminateJobObject and bounded wait on every exit, then close job."
    )
    api = _WindowsAPI()
    job = api.CreateJobObjectW(None, None)
    if not job:
        raise OSError(_windows_error("CreateJobObjectW"))
    info = None
    assigned = False
    try:
        limits = _EXTENDED_LIMITS()
        limits.BasicLimitInformation.LimitFlags = 0x00002000  # JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
        if not api.SetInformationJobObject(job, 9, ctypes.byref(limits), ctypes.sizeof(limits)):
            raise OSError(_windows_error("SetInformationJobObject(KILL_ON_JOB_CLOSE)"))
        info = _create_suspended(api, argv, cwd, stdout, stderr)
        if not api.AssignProcessToJobObject(job, info.hProcess):
            raise OSError(_windows_error("AssignProcessToJobObject (child was not resumed)"))
        assigned = True
        if time.perf_counter() >= deadline:
            completed = False
        else:
            if api.ResumeThread(info.hThread) == 0xFFFFFFFF:
                raise OSError(_windows_error("ResumeThread"))
            completed = _wait_windows(api, info.hProcess, deadline)
        if not completed:
            result["status"] = "timeout"
            result["error"] = "Process exceeded timeout_seconds."
    finally:
        try:
            exited = _cleanup_windows(api, job, info, assigned, result)
            if exited and info is not None:
                _measure_windows(api, info.hProcess, result)
        finally:
            # Keep the original process handle valid through both metric queries.
            # Closing the non-inherited job is also the last-resort tree kill.
            api.CloseHandle(job)
            if info is not None:
                api.CloseHandle(info.hThread)
                api.CloseHandle(info.hProcess)
    if result["status"] != "timeout":
        result["status"] = (
            "success" if result["returncode"] == 0
            and not result["cleanup_error"] and not result["error"] else "error"
        )
