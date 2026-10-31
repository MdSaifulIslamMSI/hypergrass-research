"""Offline checks for the optional benchmark's comparison contract."""

from pathlib import Path
import tempfile
import unittest


import run as benchmark


class FastaContractTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "input.fa"

    def parse(self, data):
        self.path.write_bytes(data)
        return benchmark.parse_fasta(self.path)

    def test_order_wrapping_line_endings_and_record_separators(self):
        expected = self.parse(b">one full header\nAcGT-N\n>two\nT.T*\n")
        actual = self.parse(b">two\r\nT.\r\nT*\r\n\r\n>one full header\r\nAc\r\nGT-N\r\n\r\n")
        benchmark.check_records(expected, actual)

    def test_invalid_fasta_is_rejected(self):
        for data in (
            b"", b"ACGT\n", b">\nACGT\n", b"> one\nACGT\n", b">one\n",
            b">one\nAC GT\n", b">one\nAC\n\nGT\n", b">one\nAC\n>one\nGT\n",
            b">one a\nAC\n>one b\nGT\n", b">one\n\nACGT\n", b">one\nAC\xff\n",
        ):
            with self.subTest(data=data), self.assertRaises(ValueError):
                self.parse(data)

    def test_record_mismatches_are_rejected_without_normalization(self):
        expected = self.parse(b">one full header\nAcGT-N\n>two\nTT\n")
        for data in (
            b">one full header\nAcGT-N\n",  # Missing record.
            b">one changed header\nAcGT-N\n>two\nTT\n",
            b">one full header\nACGT-N\n>two\nTT\n",  # Case folding is forbidden.
            b">one full header\nAcGTN\n>two\nTT\n",   # No gap removal.
            b">one full header\nAcGT-A\n>two\nTT\n", # No ambiguity replacement.
            b">one full header\nAcGT-N\n>two\nTT\n>extra\nAA\n",
        ):
            with self.subTest(data=data), self.assertRaises(ValueError):
                benchmark.check_records(expected, self.parse(data))

    def test_cached_fixture_hash_failure_does_not_replace_input(self):
        original_here, original_fixtures = benchmark.HERE, benchmark.FIXTURES
        try:
            benchmark.HERE = Path(self.directory.name)
            benchmark.FIXTURES = {"input.fa": (2, benchmark.sha256(b"AA"))}
            (benchmark.HERE / "data").mkdir()
            fixture = benchmark.HERE / "data" / "input.fa"
            fixture.write_bytes(b"CC")
            with self.assertRaises(ValueError):
                benchmark.fixtures({"fixtures": {}}, 1)
            self.assertEqual(fixture.read_bytes(), b"CC")
        finally:
            benchmark.HERE, benchmark.FIXTURES = original_here, original_fixtures

    def test_extra_size_search_is_opt_in(self):
        parser = benchmark.argument_parser()
        self.assertFalse(parser.parse_args([]).optimize)
        self.assertTrue(parser.parse_args(["--optimize"]).optimize)

    def test_pinned_demo_contract(self):
        source = benchmark.HERE / "data" / "sars_20.fa"
        if not source.exists():
            self.skipTest("Optional demo fixture has not been downloaded")
        data = source.read_bytes()
        size, digest = benchmark.FIXTURES[source.name]
        self.assertEqual((len(data), benchmark.sha256(data)), (size, digest))
        self.assertEqual(len(benchmark.parse_fasta(source)), 20)


if __name__ == "__main__":
    unittest.main()
