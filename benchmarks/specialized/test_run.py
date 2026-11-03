"""Unit tests for repeat aggregation and gate evaluation.

Run with: python -m unittest benchmarks/specialized/test_run.py -v
The harness module is loaded by path so the tests do not need a package.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


SPEC = importlib.util.spec_from_file_location(
    "specialized_run", Path(__file__).with_name("run.py")
)
specialized_run = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(specialized_run)


def observation(wall, rss, size, restore_rss=None, status="success", digest="ab", ok=True,
                cgroup_peak=None, oom_kills=None):
    entry = {
        "status": status,
        "roundtrip_ok": ok,
        "artifact_bytes": size,
        "input_bases": 1000,
        "input_bytes": 1000,
        "input_sha256": "ff",
        "wall_seconds": wall,
        "cpu_seconds": wall / 2,
        "peak_rss_bytes": rss,
        "cgroup_memory_peak_bytes": cgroup_peak,
        "cgroup_oom_kills": oom_kills,
        "verification": {"archive_sha256": digest},
        "decompress": {"peak_rss_bytes": restore_rss, "wall_seconds": wall / 3},
        "error": None if status == "success" else "boom",
    }
    return entry


class SummaryTests(unittest.TestCase):
    def profile(self):
        return {
            "id": "bounded",
            "kind": "hypergrass",
            "level": "bounded",
            "options": ["--bounded"],
            "heap": "160m",
            "memory": "256m",
        }

    def test_median_wall_and_worst_rss(self):
        summary = specialized_run.summarize_profile(
            self.profile(),
            [
                observation(20.0, 240 << 20, 1000, 200 << 20, cgroup_peak=250 << 20),
                observation(10.0, 260 << 20, 1000, 210 << 20, cgroup_peak=270 << 20, oom_kills=0),
                observation(30.0, 250 << 20, 1000, 220 << 20, cgroup_peak=265 << 20),
            ],
            3,
        )
        self.assertEqual(summary["status"], "success")
        self.assertTrue(summary["deterministic"])
        self.assertEqual(summary["artifact_bytes"], 1000)
        self.assertEqual(summary["wall_seconds"], 20.0)
        self.assertEqual(summary["wall_seconds_min"], 10.0)
        self.assertEqual(summary["wall_seconds_max"], 30.0)
        self.assertEqual(summary["peak_rss_bytes"], 260 << 20)
        self.assertEqual(summary["cgroup_memory_peak_bytes"], 270 << 20)
        self.assertEqual(summary["cgroup_oom_kills"], 0)
        self.assertEqual(summary["decompress"]["peak_rss_bytes"], 220 << 20)
        self.assertEqual(len(summary["observations"]), 3)

    def test_size_disagreement_is_not_deterministic(self):
        summary = specialized_run.summarize_profile(
            self.profile(),
            [observation(10, 1, 1000), observation(11, 2, 1001)],
            2,
        )
        self.assertEqual(summary["status"], "partial")
        self.assertFalse(summary["deterministic"])
        self.assertIsNone(summary["artifact_bytes"])
        self.assertIsNone(summary["bits_per_base"])

    def test_any_failed_repeat_makes_profile_partial(self):
        summary = specialized_run.summarize_profile(
            self.profile(),
            [
                observation(10, 1, 1000),
                observation(11, 2, None, status="oom-heap", ok=False),
            ],
            2,
        )
        self.assertEqual(summary["status"], "partial")
        self.assertFalse(summary["roundtrip_ok"])
        self.assertEqual(summary["error"], "boom")

    def test_no_successful_observation_reports_failure(self):
        summary = specialized_run.summarize_profile(
            self.profile(),
            [observation(1, 1, None, status="oom-killed", ok=False)],
            1,
        )
        self.assertEqual(summary["status"], "failed")
        self.assertIsNone(summary["wall_seconds"])
        self.assertEqual(summary["attempts"], 1)
        self.assertEqual(summary["failed_observations"], 1)
        self.assertEqual(summary["attempt_wall_seconds_median"], 1)
        self.assertEqual(summary["attempt_peak_rss_bytes"], 1)

    def test_failed_profile_keeps_attempt_cost_and_oom_evidence(self):
        summary = specialized_run.summarize_profile(
            self.profile(),
            [
                observation(10, 260 << 20, None, status="oom-killed", ok=False,
                            cgroup_peak=256 << 20, oom_kills=1),
                observation(20, 262 << 20, None, status="oom-killed", ok=False,
                            cgroup_peak=256 << 20, oom_kills=1),
                observation(30, 264 << 20, None, status="oom-killed", ok=False,
                            cgroup_peak=256 << 20, oom_kills=1),
            ],
            3,
        )
        self.assertEqual(summary["status"], "failed")
        self.assertIsNone(summary["artifact_bytes"])
        self.assertEqual(summary["attempts"], 3)
        self.assertEqual(summary["failed_observations"], 3)
        self.assertEqual(summary["attempt_wall_seconds_median"], 20)
        self.assertEqual(summary["attempt_wall_seconds_max"], 30)
        self.assertEqual(summary["attempt_peak_rss_bytes"], 264 << 20)
        self.assertEqual(summary["cgroup_memory_peak_bytes"], 256 << 20)
        self.assertEqual(summary["cgroup_oom_kills"], 1)
        self.assertEqual(summary["oom_killed_observations"], 3)
        self.assertEqual(summary["error"], "boom")


class FailureReportingTests(unittest.TestCase):
    def test_failed_restore_error_is_propagated(self):
        entry = {
            "status": "success",
            "roundtrip_ok": False,
            "decompress": {
                "status": "oom-killed",
                "error": "Container cgroup recorded an OOM kill (exit -9).",
            },
        }
        self.assertEqual(
            specialized_run.roundtrip_failure(entry),
            "Container cgroup recorded an OOM kill (exit -9).",
        )

    def test_failed_restore_without_text_is_described(self):
        entry = {
            "status": "success",
            "roundtrip_ok": False,
            "decompress": {"status": "killed", "error": None},
        }
        self.assertEqual(
            specialized_run.roundtrip_failure(entry),
            "restore failed with status killed",
        )

    def test_verification_problems_are_named(self):
        entry = {
            "status": "success",
            "roundtrip_ok": False,
            "decompress": {"status": "success"},
            "verification": {
                "verdict": "fail",
                "problems": ["restored length differs", "hash differs"],
            },
        }
        self.assertEqual(
            specialized_run.roundtrip_failure(entry),
            "independent verification failed: restored length differs; hash differs",
        )

    def test_missing_restore_is_explicit(self):
        entry = {"status": "success", "roundtrip_ok": False}
        self.assertEqual(
            specialized_run.roundtrip_failure(entry),
            "restore was not attempted",
        )


class GateTests(unittest.TestCase):
    def summary(self, profile, size, wall, rss, restore=None, status="success"):
        return {
            "profile": profile,
            "status": status,
            "artifact_bytes": size,
            "wall_seconds": wall,
            "peak_rss_bytes": rss,
            "decompress": {"peak_rss_bytes": restore if restore is not None else rss},
        }

    def test_perfect_gate_candidate_passes(self):
        gate = specialized_run.evaluate_gate(
            [
                self.summary("jarvis3-l7", 1000, 10.0, 200 << 20),
                self.summary("bounded", 1010, 20.0, 250 << 20, 240 << 20),
            ],
            "jarvis3-l7", 256.0, 1.02, 2.0,
        )
        self.assertTrue(gate["passed"], gate)
        self.assertTrue(all(check["passed"] for check in gate["checks"]))

    def test_oversized_or_slow_candidate_fails(self):
        gate = specialized_run.evaluate_gate(
            [
                self.summary("jarvis3-l7", 1000, 10.0, 200 << 20),
                self.summary("bounded", 1021, 20.1, 257 << 20, 200 << 20),
            ],
            "jarvis3-l7", 256.0, 1.02, 2.0,
        )
        self.assertFalse(gate["passed"])
        failed = {
            check["check"] for check in gate["checks"] if not check["passed"]
        }
        self.assertEqual(failed, {"compress_rss", "size_ratio", "wall_ratio"})

    def test_nondeterministic_candidate_fails_even_when_metrics_fit(self):
        candidate = self.summary("bounded", 1000, 10.0, 200 << 20)
        candidate["status"] = "partial"
        gate = specialized_run.evaluate_gate(
            [self.summary("jarvis3-l7", 1000, 10.0, 200 << 20), candidate],
            "jarvis3-l7", 256.0, 1.02, 2.0,
        )
        self.assertFalse(gate["passed"])
        self.assertIn("status", {
            check["check"] for check in gate["checks"] if not check["passed"]
        })

    def test_missing_or_failed_baseline_is_explicit_error(self):
        gate = specialized_run.evaluate_gate(
            [self.summary("bounded", 1000, 10.0, 200 << 20)],
            "jarvis3-l7", 256.0, 1.02, 2.0,
        )
        self.assertFalse(gate["passed"])
        self.assertTrue(gate["error"])
        failed = self.summary("jarvis3-l7", 1000, 10.0, 200 << 20, status="partial")
        gate = specialized_run.evaluate_gate(
            [failed, self.summary("bounded", 1000, 10.0, 200 << 20)],
            "jarvis3-l7", 256.0, 1.02, 2.0,
        )
        self.assertFalse(gate["passed"])
        self.assertTrue(gate["error"])

    def test_candidate_restriction_ignores_other_comparison_profiles(self):
        gate = specialized_run.evaluate_gate(
            [
                self.summary("jarvis3-l7", 1000, 10.0, 200 << 20),
                self.summary("bounded", 1000, 10.0, 200 << 20),
                self.summary("xz-9", 5000, 99.0, 900 << 20),
            ],
            "jarvis3-l7", 256.0, 1.02, 2.0,
            candidates=["bounded"],
        )
        self.assertTrue(gate["passed"], gate)

    def test_empty_or_unknown_candidate_list_is_an_error(self):
        summaries = [self.summary("jarvis3-l7", 1000, 10.0, 200 << 20)]
        empty = specialized_run.evaluate_gate(
            summaries, "jarvis3-l7", 256.0, 1.02, 2.0, candidates=[],
        )
        self.assertTrue(empty["error"])
        unknown = specialized_run.evaluate_gate(
            summaries, "jarvis3-l7", 256.0, 1.02, 2.0,
            candidates=["not-selected"],
        )
        self.assertIn("not-selected", unknown["error"])

    def test_unknown_rss_cannot_pass(self):
        candidate = self.summary("bounded", 1000, 10.0, None)
        candidate["decompress"] = {"peak_rss_bytes": None}
        gate = specialized_run.evaluate_gate(
            [self.summary("jarvis3-l7", 1000, 10.0, 200 << 20), candidate],
            "jarvis3-l7", 256.0, 1.02, 2.0,
        )
        self.assertFalse(gate["passed"])
        failed = {
            check["check"] for check in gate["checks"] if not check["passed"]
        }
        self.assertEqual(failed, {"compress_rss", "restore_rss"})


if __name__ == "__main__":
    unittest.main()
