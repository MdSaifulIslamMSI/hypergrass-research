"""Host-only tests; every Docker invocation is mocked."""

import copy
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

if __package__:
    from . import build_image
else:
    import build_image


class BuildImageTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="corpus-build-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        (self.root / "vendor").mkdir()
        for name in build_image.CONTEXT_FILES:
            if name != "sources.lock.json":
                shutil.copyfile(Path(build_image.__file__).parent / name, self.root / name)
        self.sources = []
        for source_id, commit in build_image.SOURCE_COMMITS.items():
            content = f"test-only archive for {source_id}".encode()
            relative = f"vendor/{source_id}.tar.gz"
            (self.root / relative).write_bytes(content)
            self.sources.append({
                "id": source_id, "commit": commit,
                "url": f"https://example.org/{source_id}/archive/{commit}.tar.gz",
                "path": relative, "size_bytes": len(content),
                "sha256": hashlib.sha256(content).hexdigest(),
            })
        self.write_sources()
        self.image_id = "sha256:" + "a" * 64
        self.environment = {
            "schema_version": 1,
            "base_image": build_image.BASE_IMAGE,
            "platform": build_image.PLATFORM,
            "alpine_release": "3.24.0",
            "sources": sorted(copy.deepcopy(self.sources), key=lambda item: item["id"]),
            "cflags": dict(build_image.CFLAGS),
            "versions": {name: f"{name} test version" for name in build_image.VERSION_COMMANDS},
            "apk_packages": ["gcc-15.2.0-r0", "musl-1.2.5-r0"],
            "executables": {
                name: {"path": f"/usr/local/bin/{name}", "sha256": "b" * 64}
                for name in build_image.EXECUTABLES
            },
        }
        self.commands = []
        self.inspection = {"Id": self.image_id, "Os": "linux", "Architecture": "amd64"}
        self.lock_path = self.root / "image.lock.json"

    def write_sources(self):
        (self.root / "sources.lock.json").write_text(
            json.dumps({"schema_version": 1, "sources": self.sources}), encoding="utf-8",
        )

    def fake_docker(self, command, **kwargs):
        self.assertIsInstance(command, list)
        self.assertNotIn("shell", kwargs)
        self.assertTrue(kwargs["check"])
        self.assertGreater(kwargs["timeout"], 0)
        self.commands.append(command)
        if command[:2] == ["docker", "build"]:
            Path(command[command.index("--iidfile") + 1]).write_text(self.image_id, encoding="utf-8")
            kwargs["stdout"].write("mock build output\n")
            return subprocess.CompletedProcess(command, 0)
        if command[:3] == ["docker", "image", "inspect"]:
            return subprocess.CompletedProcess(command, 0, json.dumps([self.inspection]))
        if command[:2] == ["docker", "run"]:
            return subprocess.CompletedProcess(command, 0, json.dumps(self.environment))
        if command[:3] == ["docker", "image", "tag"]:
            return subprocess.CompletedProcess(command, 0)
        self.fail(f"Unexpected external command: {command!r}")

    def build(self):
        with mock.patch.object(build_image.subprocess, "run", side_effect=self.fake_docker):
            with mock.patch("sys.stdout", new=io.StringIO()):
                return build_image.build_image(self.root)

    def assert_not_tagged(self):
        self.assertFalse(any(command[:3] == ["docker", "image", "tag"] for command in self.commands))

    def test_sources_validate_and_sort_without_subprocesses(self):
        with mock.patch.object(build_image.subprocess, "run") as run:
            sources = build_image.validate_sources(self.root)
        self.assertEqual([item["id"] for item in sources], ["geco3", "jarvis3"])
        run.assert_not_called()

    def test_rejects_invalid_source_fields(self):
        cases = [
            ("commit", "0" * 40, "pinned commit"),
            ("path", "../jarvis3.tar.gz", "path must be"),
            ("path", str(self.root / "vendor/jarvis3.tar.gz"), "path must be"),
            ("url", "http://example.org/source", "HTTPS"),
            ("url", None, "HTTPS"),
            ("size_bytes", True, "size_bytes"),
            ("size_bytes", -1, "size_bytes"),
            ("size_bytes", 1, "size mismatch"),
            ("sha256", "not-a-sha256", "Invalid sha256"),
            ("sha256", "0" * 64, "SHA-256 mismatch"),
        ]
        original = copy.deepcopy(self.sources)
        for field, value, message in cases:
            with self.subTest(field=field, value=value):
                self.sources = copy.deepcopy(original)
                self.sources[0][field] = value
                self.write_sources()
                with self.assertRaisesRegex(build_image.BuildError, message):
                    build_image.validate_sources(self.root)

    def test_rejects_missing_duplicate_unknown_and_malformed_sources(self):
        original = copy.deepcopy(self.sources)
        manifests = [
            {}, {"sources": None}, {"sources": []},
            {"sources": [original[0], original[0]]},
            {"sources": [dict(original[0], id="other"), original[1]]},
            {"sources": [None, original[1]]},
            {"sources": [dict(original[0], id=[]), original[1]]},
        ]
        for manifest in manifests:
            with self.subTest(manifest=manifest):
                (self.root / "sources.lock.json").write_text(json.dumps(manifest), encoding="utf-8")
                with self.assertRaises(build_image.BuildError):
                    build_image.validate_sources(self.root)

    def test_missing_or_modified_archive_rejected_before_docker(self):
        archive = self.root / self.sources[0]["path"]
        archive.unlink()
        with mock.patch.object(build_image.subprocess, "run") as run:
            with self.assertRaisesRegex(build_image.BuildError, "Missing or unsafe"):
                build_image.build_image(self.root)
            run.assert_not_called()
        archive.write_bytes(b"x" * self.sources[0]["size_bytes"])
        with mock.patch.object(build_image.subprocess, "run") as run:
            with self.assertRaisesRegex(build_image.BuildError, "SHA-256 mismatch"):
                build_image.build_image(self.root)
            run.assert_not_called()

    def test_missing_or_invalid_json_manifest_rejected(self):
        path = self.root / "sources.lock.json"
        for content in ("{", "null", "[]"):
            with self.subTest(content=content):
                path.write_text(content, encoding="utf-8")
                with self.assertRaises(build_image.BuildError):
                    build_image.validate_sources(self.root)
        path.unlink()
        with self.assertRaisesRegex(build_image.BuildError, "Cannot read"):
            build_image.validate_sources(self.root)

    def test_build_records_image_environment_and_log(self):
        result = self.build()
        self.assertEqual(json.loads(self.lock_path.read_text(encoding="utf-8")), result)
        self.assertEqual(result["image_id"], self.image_id)
        self.assertEqual(result["environment"], self.environment)
        self.assertEqual(result["build_inputs"]["base_image"], build_image.BASE_IMAGE)
        self.assertEqual(result["build_inputs"]["cflags"], build_image.CFLAGS)
        self.assertEqual(set(result["build_inputs"]["context_sha256"]), set(build_image.CONTEXT_FILES))
        log = (self.root / "build.log").read_text(encoding="utf-8")
        self.assertIn("mock build output", log)
        self.assertIn("docker build", log)
        self.assertIn("docker image tag", log)
        build_command = self.commands[0]
        self.assertEqual(build_command[-1], str(self.root.resolve()))
        self.assertEqual(build_command[build_command.index("--platform") + 1], "linux/amd64")
        self.assertNotIn("--tag", build_command)
        self.assertNotIn("-t", build_command)
        run_command = next(command for command in self.commands if command[:2] == ["docker", "run"])
        for flag in ("--rm", "--pull=never", "--network=none"):
            self.assertIn(flag, run_command)
        for flag in ("-v", "--volume", "--mount", "--privileged"):
            self.assertNotIn(flag, run_command)
        self.assertEqual(run_command[-2:], [self.image_id, "/opt/benchmark-environment.json"])
        self.assertEqual(self.commands[-1], ["docker", "image", "tag", self.image_id, build_image.IMAGE_TAG])

    def test_identical_rebuild_preserves_lock_bytes_and_appends_log(self):
        self.build()
        original = self.lock_path.read_bytes()
        original_size = (self.root / "build.log").stat().st_size
        self.build()
        self.assertEqual(self.lock_path.read_bytes(), original)
        self.assertGreater((self.root / "build.log").stat().st_size, original_size)

    def test_changed_context_refuses_existing_lock_before_docker(self):
        self.build()
        original = self.lock_path.read_bytes()
        with (self.root / "Dockerfile").open("a", encoding="utf-8") as stream:
            stream.write("\n# changed build recipe\n")
        self.commands.clear()
        with self.assertRaisesRegex(build_image.BuildError, "Incompatible image.lock"):
            self.build()
        self.assertEqual(self.commands, [])
        self.assertEqual(self.lock_path.read_bytes(), original)

    def test_malformed_or_incompatible_lock_is_not_overwritten(self):
        for content in ("{", "null", "[]", "{}", '{"schema_version": 99}'):
            with self.subTest(content=content):
                self.lock_path.write_text(content, encoding="utf-8")
                with self.assertRaises(build_image.BuildError):
                    self.build()
                self.assertEqual(self.lock_path.read_text(encoding="utf-8"), content)
                self.assertEqual(self.commands, [])

    def test_changed_image_id_refuses_to_retag_or_overwrite_lock(self):
        self.build()
        original = self.lock_path.read_bytes()
        self.commands.clear()
        self.image_id = "sha256:" + "c" * 64
        self.inspection["Id"] = self.image_id
        with self.assertRaisesRegex(build_image.BuildError, "Rebuilt image differs"):
            self.build()
        self.assertEqual(self.lock_path.read_bytes(), original)
        self.assert_not_tagged()

    def test_changed_package_versions_refuses_to_retag_or_overwrite_lock(self):
        self.build()
        original = self.lock_path.read_bytes()
        self.commands.clear()
        self.environment["apk_packages"].append("new-package-1.0-r0")
        with self.assertRaisesRegex(build_image.BuildError, "Rebuilt image differs"):
            self.build()
        self.assertEqual(self.lock_path.read_bytes(), original)
        self.assert_not_tagged()

    def test_wrong_platform_or_image_id_rejected(self):
        for field, value in (("Architecture", "arm64"), ("Os", "windows"), ("Id", "sha256:" + "d" * 64)):
            with self.subTest(field=field):
                self.inspection = {"Id": self.image_id, "Os": "linux", "Architecture": "amd64"}
                self.inspection[field] = value
                self.commands.clear()
                with self.assertRaisesRegex(build_image.BuildError, "linux/amd64"):
                    self.build()
                self.assertFalse(self.lock_path.exists())
                self.assert_not_tagged()

    def test_environment_must_include_correct_pins_versions_and_hashes(self):
        invalid = [
            ("cflags", {"jarvis3": "-march=native"}),
            ("base_image", "eclipse-temurin:latest"),
            ("sources", []), ("alpine_release", "3.23.0"),
            ("versions", {}), ("apk_packages", []), ("executables", {}),
        ]
        for key, value in invalid:
            with self.subTest(key=key):
                environment = dict(self.environment, **{key: value})
                with self.assertRaises(build_image.BuildError):
                    build_image.validate_environment(environment, self.environment["sources"])
        for field, value in (("sha256", "bad"), ("path", "/tmp/JARVIS3")):
            with self.subTest(field=field):
                environment = copy.deepcopy(self.environment)
                environment["executables"]["JARVIS3"][field] = value
                with self.assertRaises(build_image.BuildError):
                    build_image.validate_environment(environment, self.environment["sources"])

    def test_bad_environment_never_creates_lock_or_tag(self):
        self.environment["cflags"] = {}
        with self.assertRaisesRegex(build_image.BuildError, "incompatible cflags"):
            self.build()
        self.assertFalse(self.lock_path.exists())
        self.assert_not_tagged()

    def test_build_failure_keeps_output_without_creating_lock(self):
        def fail(command, **kwargs):
            kwargs["stdout"].write("compiler failure\n")
            raise subprocess.CalledProcessError(2, command)

        with mock.patch.object(build_image.subprocess, "run", side_effect=fail):
            with mock.patch("sys.stdout", new=io.StringIO()):
                with self.assertRaises(subprocess.CalledProcessError):
                    build_image.build_image(self.root)
        self.assertIn("compiler failure", (self.root / "build.log").read_text(encoding="utf-8"))
        self.assertFalse(self.lock_path.exists())

    def test_inputs_changed_during_build_are_not_locked_or_tagged(self):
        def change_input(command, **kwargs):
            result = self.fake_docker(command, **kwargs)
            if command[:2] == ["docker", "build"]:
                with (self.root / "Dockerfile").open("a", encoding="utf-8") as stream:
                    stream.write("\n# concurrent edit\n")
            return result

        with mock.patch.object(build_image.subprocess, "run", side_effect=change_input):
            with mock.patch("sys.stdout", new=io.StringIO()):
                with self.assertRaisesRegex(build_image.BuildError, "inputs changed"):
                    build_image.build_image(self.root)
        self.assertFalse(self.lock_path.exists())
        self.assert_not_tagged()

    def test_environment_recorder_captures_versions_packages_and_executable_hashes(self):
        output = self.root / "environment.json"
        executable = self.root / "test-executable"
        executable.write_bytes(b"test executable")
        original_read_text = Path.read_text

        def read_text(path, *args, **kwargs):
            if path == Path("/etc/alpine-release"):
                return "3.24.0\n"
            return original_read_text(path, *args, **kwargs)

        def run(command, **kwargs):
            self.assertIsInstance(command, list)
            if command == ["apk", "info", "-v"]:
                stdout = "musl-1.2.5-r0\ngcc-15.2.0-r0\n"
            else:
                stdout = "test version output\n"
                self.assertEqual(kwargs["stderr"], subprocess.STDOUT)
            return subprocess.CompletedProcess(command, 0, stdout)

        digest = hashlib.sha256(executable.read_bytes()).hexdigest()
        with (
            mock.patch.object(build_image, "validate_sources", return_value=self.environment["sources"]),
            mock.patch.dict(build_image.os.environ, {
                "JARVIS3_CFLAGS": build_image.CFLAGS["jarvis3"],
                "GECO3_CFLAGS": build_image.CFLAGS["geco3"],
            }),
            mock.patch.object(build_image.shutil, "which", side_effect=lambda name: f"/usr/local/bin/{name}"),
            mock.patch.object(build_image, "sha256_file", return_value=digest),
            mock.patch.object(Path, "read_text", new=read_text),
            mock.patch.object(build_image.subprocess, "run", side_effect=run) as run_mock,
        ):
            build_image.record_environment(output)
        environment = json.loads(output.read_text(encoding="utf-8"))
        self.assertEqual(environment["versions"]["java"], "test version output")
        self.assertEqual(environment["apk_packages"], ["gcc-15.2.0-r0", "musl-1.2.5-r0"])
        self.assertEqual(environment["executables"]["JARVIS3"]["sha256"], hashlib.sha256(executable.read_bytes()).hexdigest())
        self.assertEqual(run_mock.call_count, len(build_image.VERSION_COMMANDS) + 1)

    def test_internal_verify_mode_never_calls_docker(self):
        with mock.patch.object(build_image.subprocess, "run") as run:
            self.assertEqual(build_image.main(["--verify-sources", str(self.root)]), 0)
        run.assert_not_called()

    def test_dockerfile_and_context_are_pinned_and_minimal(self):
        dockerfile = (self.root / "Dockerfile").read_text(encoding="utf-8")
        self.assertTrue(dockerfile.startswith(f"FROM --platform=linux/amd64 {build_image.BASE_IMAGE}\n"))
        self.assertIn("apk add --no-cache python3 build-base gzip xz bash", dockerfile)
        self.assertIn('JARVIS3_CFLAGS="-Wall -Ofast -DPROGRESS"', dockerfile)
        self.assertIn('GECO3_CFLAGS="-Wall -O3"', dockerfile)
        self.assertIn('make -C /opt/sources/jarvis3/src CC=gcc CFLAGS="$JARVIS3_CFLAGS"', dockerfile)
        self.assertIn('make -C /opt/sources/geco3/src CC=gcc CFLAGS="$GECO3_CFLAGS"', dockerfile)
        self.assertNotIn("-march=native", dockerfile)
        self.assertIn("RUN --network=none", dockerfile)
        self.assertIn("WORKDIR /work", dockerfile)
        self.assertIn("--record-environment /opt/benchmark-environment.json", dockerfile)
        ignores = (self.root / ".dockerignore").read_text(encoding="utf-8").splitlines()
        self.assertEqual(ignores, [
            "**", "!Dockerfile", "!.dockerignore", "!build_image.py", "!sources.lock.json",
            "!vendor/", "vendor/**", "!vendor/jarvis3.tar.gz", "!vendor/geco3.tar.gz",
        ])


if __name__ == "__main__":
    unittest.main()
