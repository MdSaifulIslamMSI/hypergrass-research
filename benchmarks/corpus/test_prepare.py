"""Offline stdlib tests: python -B -m unittest discover -s benchmarks/corpus -p test_prepare.py -v."""

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock
from urllib.parse import parse_qs, urlsplit
from urllib.request import Request

if __package__:
    from . import prepare
else:
    import prepare


A = "NC_000001.1"
B = "NC_000002.2"
FASTA_A = b">NC_000001.1 First full header  \r\nACgtNRY-.\r\n ?ac\r\n"
FASTA_B = b">NC_000002.2 Second full header\nACGT\nTGCA\n"
SEQUENCE_A = b"ACgtNRY-. ?ac"
SEQUENCE_B = b"ACGTTGCA"


class PrepareTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="corpus-prepare-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.protocol = {
            "schema_version": 1,
            "datasets": [
                {"id": "one", "label": "Mixed case and ambiguity", "accessions": [A]},
                {"id": "two", "label": "Uppercase ACGT", "accessions": [B]},
                {"id": "pair", "label": "Reverse order", "accessions": [B, A]},
            ],
        }
        self.write_protocol()

    def write_protocol(self):
        (self.root / "protocol.json").write_text(json.dumps(self.protocol), encoding="ascii")

    def download(self, url, max_bytes):
        self.assertEqual(urlsplit(url).scheme, "https")
        query = parse_qs(urlsplit(url).query)
        self.assertEqual(query["db"], ["nuccore"])
        self.assertEqual(query["rettype"], ["fasta"])
        self.assertEqual(query["retmode"], ["text"])
        self.assertEqual(max_bytes, prepare.FASTA_LIMIT)
        return {A: FASTA_A, B: FASTA_B}[query["id"][0]]

    def prepare_data(self):
        with mock.patch.object(prepare, "_download", side_effect=self.download) as download:
            manifest = prepare.prepare_data(self.root)
        self.assertEqual(download.call_count, 2)
        self.assertEqual(
            [parse_qs(urlsplit(call.args[0]).query)["id"][0] for call in download.call_args_list],
            [A, B],
        )
        return manifest

    def test_exact_bytes_headers_case_ambiguity_and_requested_order(self):
        manifest = self.prepare_data()
        self.assertEqual(set(manifest), {"schema_version", "protocol_sha256", "datasets"})
        self.assertEqual(manifest["schema_version"], 1)
        self.assertEqual(manifest["protocol_sha256"], prepare.sha256_file(self.root / "protocol.json"))
        expected = [
            (FASTA_A, SEQUENCE_A, [A]),
            (FASTA_B, SEQUENCE_B, [B]),
            (FASTA_B + FASTA_A, SEQUENCE_B + SEQUENCE_A, [B, A]),
        ]
        for entry, (fasta, sequence, accessions) in zip(manifest["datasets"], expected):
            with self.subTest(dataset=entry["id"]):
                self.assertEqual(set(entry), {"id", "label", "accessions", "fasta", "sequence", "records"})
                self.assertEqual(entry["accessions"], accessions)
                self.assertEqual([record["id"] for record in entry["records"]], accessions)
                for name, data in (("fasta", fasta), ("sequence", sequence)):
                    descriptor = entry[name]
                    self.assertEqual((self.root / descriptor["path"]).read_bytes(), data)
                    self.assertEqual(descriptor["size_bytes"], len(data))
                    self.assertEqual(descriptor["sha256"], prepare.sha256_bytes(data))
                self.assertEqual(entry["sequence"]["raw_acgt"], entry["id"] == "two")
        record = manifest["datasets"][0]["records"][0]
        self.assertEqual(record, {
            "id": A, "header": "NC_000001.1 First full header  ", "sequence_bytes": len(SEQUENCE_A),
        })
        self.assertEqual(prepare.verify_data(self.root), manifest)

    def test_sha256_helpers(self):
        path = self.root / "hash.bin"
        path.write_bytes(b"abc")
        expected = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        self.assertEqual(prepare.sha256_bytes(b"abc"), expected)
        self.assertEqual(prepare.sha256_file(path), expected)

    def test_parser_line_endings_and_no_final_newline(self):
        sequence, records = prepare._parse_fasta(b">NC_000001.1 title\rAcN\r\n-?.\nacgt", [A])
        self.assertEqual(sequence, b"AcN-?.acgt")
        self.assertEqual(records[0]["sequence_bytes"], 10)

    def test_missing_extra_wrong_version_duplicate_and_out_of_order_records(self):
        cases = [
            (b"", [A]),
            (FASTA_A, [A, B]),
            (FASTA_A + FASTA_B, [A]),
            (FASTA_B + FASTA_A, [A, B]),
            (FASTA_A.replace(b"NC_000001.1", b"NC_000001.2"), [A]),
            (FASTA_A + FASTA_A, [A, B]),
            (FASTA_A + FASTA_A.replace(b"First full header", b"Different header"), [A, B]),
        ]
        for fasta, accessions in cases:
            with self.subTest(fasta=fasta, accessions=accessions):
                with self.assertRaisesRegex(ValueError, "IDs/order|Duplicate"):
                    prepare._parse_fasta(fasta, accessions)

    def test_nonprintable_sequence_and_invalid_headers_are_rejected_not_stripped(self):
        for invalid in (b"\x00", b"\t", b"\x0b", b"\x0c", b"\x1c", b"\x7f", b"\x80"):
            with self.subTest(byte=invalid):
                with self.assertRaisesRegex(ValueError, "non-printable"):
                    prepare._parse_fasta(b">NC_000001.1 title\nAC" + invalid + b"GT\n", [A])
        for fasta in (b"ACGT\n" + FASTA_A, b">\nACGT\n", b"> NC_000001.1\nACGT\n",
                      b">NC_000001.1\tbad\nACGT\n", b">NC_000001.1\n"):
            with self.subTest(fasta=fasta):
                with self.assertRaises(ValueError):
                    prepare._parse_fasta(fasta, [A])

    def test_individual_response_must_be_exact_requested_record(self):
        for response in (FASTA_B, FASTA_A + FASTA_B, FASTA_A + FASTA_A):
            with self.subTest(response=response):
                with mock.patch.object(prepare, "_download", return_value=response) as download:
                    with self.assertRaises(ValueError):
                        prepare.prepare_data(self.root)
                self.assertEqual(download.call_count, 1)
                self.assertFalse((self.root / "data.lock.json").exists())
                self.assertFalse((self.root / "data").exists())

    def test_unseparated_records_cannot_be_silently_rewrapped(self):
        def download(url, limit):
            fasta = self.download(url, limit)
            return fasta.rstrip(b"\r\n") if fasta == FASTA_B else fasta
        with mock.patch.object(prepare, "_download", side_effect=download):
            with self.assertRaisesRegex(ValueError, "IDs/order"):
                prepare.prepare_data(self.root)
        self.assertFalse((self.root / "data.lock.json").exists())

    def test_cached_data_and_lock_are_unchanged_without_network(self):
        manifest = self.prepare_data()
        paths = list((self.root / "data").iterdir()) + [self.root / "data.lock.json"]
        before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in paths}
        with mock.patch.object(prepare, "_download", side_effect=AssertionError("network")):
            self.assertEqual(prepare.prepare_data(self.root), manifest)
            self.assertEqual(prepare.verify_data(self.root), manifest)
        self.assertEqual(before, {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in paths})

    def test_corrupt_cached_fasta_size_hash_and_sequence_rejected_without_repair(self):
        self.prepare_data()
        lock_bytes = (self.root / "data.lock.json").read_bytes()
        for relative, change in (
            ("data/one.fasta", lambda data: data.replace(b"ACgt", b"TCgt")),
            ("data/one.fasta", lambda data: data.replace(b"First", b"Other")),
            ("data/one.fasta", lambda data: data + b"\n"),
            ("data/one.seq", lambda data: data.replace(b"ACgt", b"TCgt")),
            ("data/one.seq", lambda data: data + b"A"),
        ):
            path = self.root / relative
            original = path.read_bytes()
            corrupt = change(original)
            with self.subTest(path=relative, data=corrupt):
                path.write_bytes(corrupt)
                with mock.patch.object(prepare, "_download", side_effect=AssertionError("network")):
                    with self.assertRaises(ValueError):
                        prepare.verify_data(self.root)
                    with self.assertRaises(ValueError):
                        prepare.prepare_data(self.root)
                self.assertEqual(path.read_bytes(), corrupt)
                self.assertEqual((self.root / "data.lock.json").read_bytes(), lock_bytes)
                path.write_bytes(original)

    def test_missing_cached_file_is_not_redownloaded(self):
        self.prepare_data()
        (self.root / "data/one.fasta").unlink()
        with mock.patch.object(prepare, "_download", side_effect=AssertionError("network")):
            with self.assertRaises(FileNotFoundError):
                prepare.prepare_data(self.root)

    def test_protocol_byte_change_is_rejected_even_if_semantically_identical(self):
        self.prepare_data()
        lock = (self.root / "data.lock.json").read_bytes()
        with (self.root / "protocol.json").open("ab") as stream:
            stream.write(b"\n")
        with mock.patch.object(prepare, "_download", side_effect=AssertionError("network")):
            with self.assertRaisesRegex(ValueError, "protocol_sha256"):
                prepare.prepare_data(self.root)
        self.assertEqual((self.root / "data.lock.json").read_bytes(), lock)

    def test_changed_manifest_metadata_is_rejected_and_not_replaced(self):
        manifest = self.prepare_data()
        lock_path = self.root / "data.lock.json"
        for change in (
            lambda value: value.update(schema_version=2),
            lambda value: value["datasets"][0]["fasta"].update(size_bytes=1),
            lambda value: value["datasets"][0]["fasta"].update(sha256="0" * 64),
            lambda value: value["datasets"][0]["fasta"].update(path="../outside"),
            lambda value: value["datasets"][0]["records"][0].update(header="changed"),
            lambda value: value["datasets"][0]["records"][0].update(sequence_bytes=0),
            lambda value: value["datasets"][0]["sequence"].update(raw_acgt=True),
            lambda value: value["datasets"].reverse(),
        ):
            changed = json.loads(json.dumps(manifest))
            change(changed)
            lock_path.write_text(json.dumps(changed), encoding="ascii")
            before = lock_path.read_bytes()
            with self.subTest(changed=changed):
                with self.assertRaisesRegex(ValueError, "Immutable"):
                    prepare.prepare_data(self.root)
                self.assertEqual(lock_path.read_bytes(), before)

    def test_unlocked_existing_artifacts_are_not_overwritten(self):
        data_dir = self.root / "data"
        data_dir.mkdir()
        path = data_dir / "one.fasta"
        path.write_bytes(b"untrusted old bytes")
        with mock.patch.object(prepare, "_download", side_effect=self.download):
            with self.assertRaisesRegex(ValueError, "Refusing to replace"):
                prepare.prepare_data(self.root)
        self.assertEqual(path.read_bytes(), b"untrusted old bytes")
        self.assertFalse((self.root / "data.lock.json").exists())

    def test_protocol_rejects_unpinned_and_duplicate_accessions_and_unsafe_ids(self):
        invalid = [
            {"id": "one", "label": "bad", "accessions": ["NC_000001"]},
            {"id": "one", "label": "bad", "accessions": [A, A]},
            {"id": "../escape", "label": "bad", "accessions": [A]},
        ]
        for spec in invalid:
            self.protocol["datasets"] = [spec]
            self.write_protocol()
            with self.subTest(spec=spec):
                with mock.patch.object(prepare, "_download", side_effect=AssertionError("network")):
                    with self.assertRaises(ValueError):
                        prepare.prepare_data(self.root)

    def test_sources_pin_urls_schema_cache_and_corruption(self):
        with mock.patch.object(prepare, "_download", side_effect=[b"jarvis archive", b"geco archive"]) as download:
            manifest = prepare.prepare_sources(self.root)
        self.assertEqual(set(manifest), {"schema_version", "sources"})
        self.assertEqual(manifest["schema_version"], 1)
        expected_sources = (
            ("jarvis3", "4dcb37b05b1a1482c2bf864d4968244817b0ad0a"),
            ("geco3", "69642a0a913e543599f35aee372e300601e4cbbe"),
        )
        for entry, (source_id, commit), call in zip(manifest["sources"], expected_sources, download.call_args_list):
            self.assertEqual(set(entry), {"id", "commit", "url", "path", "size_bytes", "sha256"})
            self.assertEqual(entry["id"], source_id)
            self.assertEqual(entry["commit"], commit)
            self.assertEqual(entry["url"], f"https://codeload.github.com/cobilab/{source_id}/tar.gz/{commit}")
            self.assertEqual(call.args, (entry["url"], prepare.SOURCE_LIMIT))
            self.assertEqual(entry["path"], f"vendor/{source_id}.tar.gz")
            self.assertEqual(entry["sha256"], prepare.sha256_file(self.root / entry["path"]))
        lock = self.root / "sources.lock.json"
        before = lock.read_bytes(), lock.stat().st_mtime_ns
        with mock.patch.object(prepare, "_download", side_effect=AssertionError("network")):
            self.assertEqual(prepare.prepare_sources(self.root), manifest)
            path = self.root / "vendor/jarvis3.tar.gz"
            original = path.read_bytes()
            for corrupt in (original + b"!", b"X" + original[1:]):
                path.write_bytes(corrupt)
                with self.assertRaisesRegex(ValueError, "Immutable"):
                    prepare.prepare_sources(self.root)
                self.assertEqual(path.read_bytes(), corrupt)
            path.write_bytes(original)
            changed = json.loads(lock.read_bytes())
            changed["sources"][0]["commit"] = "0" * 40
            lock.write_text(json.dumps(changed), encoding="ascii")
            with self.assertRaisesRegex(ValueError, "Immutable"):
                prepare.prepare_sources(self.root)
            lock.write_bytes(before[0])
            path.unlink()
            with self.assertRaises(FileNotFoundError):
                prepare.prepare_sources(self.root)
        self.assertEqual(lock.read_bytes(), before[0])

    def test_cli_default_modes_and_offline_verification(self):
        fake_file = self.root / "prepare.py"
        for argv, data_calls, source_calls, verify_calls in (
            ([], 1, 1, 0), (["--data-only"], 1, 0, 0),
            (["--sources-only"], 0, 1, 0), (["--verify"], 0, 0, 1),
        ):
            with self.subTest(argv=argv), contextlib.ExitStack() as stack:
                stack.enter_context(mock.patch.object(prepare, "__file__", str(fake_file)))
                stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
                data = stack.enter_context(mock.patch.object(prepare, "prepare_data", return_value={"datasets": []}))
                sources = stack.enter_context(mock.patch.object(prepare, "prepare_sources", return_value={"sources": []}))
                verify = stack.enter_context(mock.patch.object(prepare, "verify_data", return_value={"datasets": []}))
                verify_sources = stack.enter_context(mock.patch.object(prepare, "_verify_sources", return_value={"sources": []}))
                self.assertEqual(prepare.main(argv), 0)
                self.assertEqual(data.call_count, data_calls)
                self.assertEqual(sources.call_count, source_calls)
                self.assertEqual(verify.call_count, verify_calls)
                self.assertEqual(verify_sources.call_count, verify_calls)
                for function in (data, sources, verify, verify_sources):
                    if function.called:
                        function.assert_called_once_with(self.root)
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as error:
                prepare.main(["--sources-only", "--data-only"])
        self.assertEqual(error.exception.code, 2)


class DownloadTests(unittest.TestCase):
    def response(self, body, headers=None):
        response = mock.MagicMock()
        response.__enter__.return_value = response
        response.geturl.return_value = "https://example.org/file"
        response.headers = headers or {}
        stream = io.BytesIO(body)
        response.read1.side_effect = stream.read1
        return response

    def test_timeout_rate_limit_and_exact_bounded_reads(self):
        response = self.response(FASTA_A, {"Content-Length": str(len(FASTA_A))})
        opener = mock.Mock()
        opener.open.return_value = response
        with mock.patch.object(prepare, "_last_request", 100.0), \
                mock.patch.object(prepare.time, "monotonic", return_value=100.1), \
                mock.patch.object(prepare.time, "sleep") as sleep, \
                mock.patch.object(prepare, "build_opener", return_value=opener):
            self.assertEqual(prepare._download("https://example.org/file", len(FASTA_A)), FASTA_A)
        sleep.assert_called_once()
        self.assertAlmostEqual(sleep.call_args.args[0], prepare.REQUEST_INTERVAL - 0.1)
        self.assertEqual(opener.open.call_args.kwargs, {"timeout": 30})
        request = opener.open.call_args.args[0]
        self.assertEqual(request.get_header("Accept-encoding"), "identity")
        for call in response.read1.call_args_list:
            self.assertGreater(call.args[0], 0)
            self.assertLessEqual(call.args[0], prepare.READ_CHUNK)

    def test_oversized_truncated_encoded_and_insecure_responses_are_rejected(self):
        cases = [
            (b"123456", {}, "https://example.org/file"),
            (b"12", {"Content-Length": "6"}, "https://example.org/file"),
            (b"12", {"Content-Length": "3"}, "https://example.org/file"),
            (b"12", {"Content-Encoding": "gzip"}, "https://example.org/file"),
            (b"12", {}, "http://example.org/file"),
        ]
        for body, headers, final_url in cases:
            response = self.response(body, headers)
            response.geturl.return_value = final_url
            opener = mock.Mock()
            opener.open.return_value = response
            with self.subTest(body=body, headers=headers, final_url=final_url), \
                    mock.patch.object(prepare, "_last_request", None), \
                    mock.patch.object(prepare, "build_opener", return_value=opener):
                with self.assertRaises(ValueError):
                    prepare._download("https://example.org/file", 5)
        with self.assertRaisesRegex(ValueError, "HTTPS"):
            prepare._download("http://example.org/file", 5)
        handler = prepare._HTTPSRedirectHandler()
        with self.assertRaisesRegex(ValueError, "non-HTTPS"):
            handler.redirect_request(Request("https://example.org/file"), None, 302, "Found", {}, "http://example.org/file")

    def test_elapsed_download_deadline(self):
        opener = mock.Mock()
        opener.open.return_value = self.response(b"abc")
        with mock.patch.object(prepare, "_last_request", None), \
                mock.patch.object(prepare.time, "monotonic", side_effect=[0.0, 30.1]), \
                mock.patch.object(prepare, "build_opener", return_value=opener):
            with self.assertRaises(TimeoutError):
                prepare._download("https://example.org/file", 5)


if __name__ == "__main__":
    unittest.main()
