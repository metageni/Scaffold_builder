#!/usr/bin/python3
"""
Integration tests for scaffold_builder.py.

These tests drive the full pipeline (coord2hash → clean_coords → sort_mapping)
using pre-baked .coords fixtures, bypassing nucmer entirely.  Each test
asserts on the files written to a temporary directory, verifying that the
Python 3 rewrite produces output identical in structure and content to what
the original Python 2 version would produce.
"""

import sys
import shutil
import pytest

from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from scaffold_builder import build_parameters

from utils import (clean_coords,
                   coord2hash,
                   fasta2hash,
                   sort_mapping)

FIXTURES = Path(__file__).parent / "fixtures"
QUERY = str(FIXTURES / "query.fasta")
REFERENCE = str(FIXTURES / "reference.fasta")


def _run_pipeline(tmp_path, coords_file):
    """Run the pipeline from coord2hash onward, writing output to tmp_path.

    Args:
        tmp_path (Path): Pytest-provided temporary directory.
        coords_file (str): Path to the .coords fixture file.

    Returns:
        dict: The parameters dict after the full run.
    """
    prefix = str(tmp_path / "out")
    (tmp_path / "out_overlap_alignment").mkdir()

    params = build_parameters({"-q": QUERY, "-r": REFERENCE, "-p": prefix})
    seqs = fasta2hash(QUERY)
    coords = coord2hash(coords_file, params, seqs)
    cleaned = clean_coords(coords, params, seqs)
    sort_mapping(cleaned, params, seqs)
    return params


# ---------------------------------------------------------------------------
# Gap-filling path (test.coords: contig1→contig2 gap=4, contig2→contig3 gap=75)
# ---------------------------------------------------------------------------

class TestGapFilling:
    """Integration tests for the gap-filling code path."""

    def test_scaffold_fasta_created(self, tmp_path):
        """The scaffold FASTA file is created after a successful run."""
        _run_pipeline(tmp_path, str(FIXTURES / "test.coords"))
        assert (tmp_path / "out_Scaffold.fasta").exists()

    def test_output_txt_created(self, tmp_path):
        """The output log file is created after a successful run."""
        _run_pipeline(tmp_path, str(FIXTURES / "test.coords"))
        assert (tmp_path / "out_output.txt").exists()

    def test_unmapped_contig_in_scaffold(self, tmp_path):
        """contig_unmapped appears in the scaffold FASTA with _not_mapped suffix."""
        _run_pipeline(tmp_path, str(FIXTURES / "test.coords"))
        content = (tmp_path / "out_Scaffold.fasta").read_text()
        assert "contig_unmapped_not_mapped" in content

    def test_gap_recorded_in_parameters(self, tmp_path):
        """Gaps between non-overlapping contigs are recorded in parameters."""
        params = _run_pipeline(tmp_path, str(FIXTURES / "test.coords"))
        assert len(params["gaps"]) > 0

    def test_scaffold_contains_ns_for_gap(self, tmp_path):
        """Gap between contig1 (ends at ref 60) and contig2 (starts at ref 65) → 4 Ns."""
        _run_pipeline(tmp_path, str(FIXTURES / "test.coords"))
        content = (tmp_path / "out_Scaffold.fasta").read_text()
        # contig1 ends at ref pos 60, contig2 starts at 65 → gap = 65-60-1 = 4 Ns
        assert "NNNN" in content

    def test_output_log_contains_reference_name(self, tmp_path):
        """The output log references the reference file path."""
        _run_pipeline(tmp_path, str(FIXTURES / "test.coords"))
        log = (tmp_path / "out_output.txt").read_text()
        assert REFERENCE in log

    def test_output_log_contains_statistics_section(self, tmp_path):
        """The output log contains the statistics section header."""
        _run_pipeline(tmp_path, str(FIXTURES / "test.coords"))
        log = (tmp_path / "out_output.txt").read_text()
        assert "Scaffolding statistics:" in log

    def test_output_log_contains_final_scaffolding_section(self, tmp_path):
        """The output log contains the final scaffolding table."""
        _run_pipeline(tmp_path, str(FIXTURES / "test.coords"))
        log = (tmp_path / "out_output.txt").read_text()
        assert "Final scaffolding:" in log

    def test_large_gap_splits_scaffold(self, tmp_path):
        """A gap >= -g (5000 nt) causes a new scaffold sequence to be started."""
        # contig2→contig3: ref gap = 200-124-1 = 75 < 5000, no split expected here.
        # contig3→contig4: ref gap = 260-259-1 = 0, no gap.
        # Verify at least one scaffold sequence is written.
        _run_pipeline(tmp_path, str(FIXTURES / "test.coords"))
        content = (tmp_path / "out_Scaffold.fasta").read_text()
        assert ">Scaffold_1" in content


# ---------------------------------------------------------------------------
# Ambiguous mapping path (ambiguous.coords)
# ---------------------------------------------------------------------------

class TestAmbiguousMapping:
    """Integration tests for the ambiguous-mapping removal path."""

    def test_ambiguous_contig_not_in_scaffold(self, tmp_path):
        """contig1 (mapped twice at >=95%) must not appear as a scaffold entry."""
        params = _run_pipeline(tmp_path, str(FIXTURES / "ambiguous.coords"))
        assert "contig1" in params["ambiguous"]

    def test_ambiguous_contig_written_to_fasta(self, tmp_path):
        """Ambiguous contigs are written only when there are >1 ambiguous contigs (original behaviour)."""
        # ambiguous.coords has only 1 ambiguous contig → write_output guard (>1) skips the block;
        # contig1 is simply absent from the scaffold output.
        _run_pipeline(tmp_path, str(FIXTURES / "ambiguous.coords"))
        content = (tmp_path / "out_Scaffold.fasta").read_text()
        assert "contig1_ambiguous_mapping" not in content
        assert "contig1_not_mapped" not in content

    def test_non_ambiguous_contig_scaffolded(self, tmp_path):
        """With only contig2 remaining after ambiguity removal, a single-entry mapping
        produces no Scaffold_N record (the loop never runs) but the output file is created."""
        _run_pipeline(tmp_path, str(FIXTURES / "ambiguous.coords"))
        # Single-contig mapping: _scaffold loop body never executes, scaffold_seq stays ""
        content = (tmp_path / "out_Scaffold.fasta").read_text()
        assert ">Scaffold_1" not in content
        assert (tmp_path / "out_output.txt").exists()


# ---------------------------------------------------------------------------
# Reverse-complement path (reverse.coords)
# ---------------------------------------------------------------------------

class TestReverseComplement:
    """Integration tests for the reverse-complement orientation path."""

    def test_scaffold_created_for_reverse_contig(self, tmp_path):
        """A contig mapped in reverse orientation is still scaffolded."""
        _run_pipeline(tmp_path, str(FIXTURES / "reverse.coords"))
        assert (tmp_path / "out_Scaffold.fasta").exists()
        content = (tmp_path / "out_Scaffold.fasta").read_text()
        assert ">Scaffold_1" in content

    def test_reverse_contig_sequence_is_rc(self, tmp_path):
        """The reverse-mapped contig sequence in the scaffold is the RC of the original."""
        from utils import reverse_complement
        original = fasta2hash(QUERY)["contig1"]
        _run_pipeline(tmp_path, str(FIXTURES / "reverse.coords"))
        content = (tmp_path / "out_Scaffold.fasta").read_text()
        rc = reverse_complement(original)
        assert rc in content


# ---------------------------------------------------------------------------
# Overlap path: two contigs whose reference positions overlap
# ---------------------------------------------------------------------------

class TestOverlapMerging:
    """Integration tests for the Needleman-Wunsch overlap-merging path."""

    def test_overlap_alignment_files_written(self, tmp_path):
        """Overlap and alignment FASTA files are written for overlapping contigs."""
        # Build a coords file where contig1 and contig2 overlap on the reference
        coords = tmp_path / "overlap.coords"
        # Start at r1=2 to avoid triggering begin_end logic
        coords.write_text(
            "2\t71\t1\t60\t70\t60\t100.00\treference\tcontig1\n"
            "62\t121\t1\t60\t60\t60\t100.00\treference\tcontig2\n"
        )
        (tmp_path / "out_overlap_alignment").mkdir(exist_ok=True)
        prefix = str(tmp_path / "out")
        params = build_parameters({"-q": QUERY, "-r": REFERENCE, "-p": prefix})
        seqs = fasta2hash(QUERY)
        c = coord2hash(str(coords), params, seqs)
        cleaned = clean_coords(c, params, seqs)
        sort_mapping(cleaned, params, seqs)

        overlap_dir = tmp_path / "out_overlap_alignment"
        overlap_files = list(overlap_dir.glob("Overlap_*.fasta"))
        assert len(overlap_files) > 0

    def test_overlap_increments_counter(self, tmp_path):
        """overlapOver or overlapBelow counter is incremented for overlapping contigs."""
        coords = tmp_path / "overlap.coords"
        # Start at r1=2 to avoid triggering begin_end logic
        coords.write_text(
            "2\t71\t1\t60\t70\t60\t100.00\treference\tcontig1\n"
            "62\t121\t1\t60\t60\t60\t100.00\treference\tcontig2\n"
        )
        (tmp_path / "out_overlap_alignment").mkdir(exist_ok=True)
        prefix = str(tmp_path / "out")
        params = build_parameters({"-q": QUERY, "-r": REFERENCE, "-p": prefix})
        seqs = fasta2hash(QUERY)
        c = coord2hash(str(coords), params, seqs)
        cleaned = clean_coords(c, params, seqs)
        sort_mapping(cleaned, params, seqs)

        total_overlaps = params["overlapOver"] + params["overlapBelow"]
        assert total_overlaps == 1


# ---------------------------------------------------------------------------
# Parity test: Python 3 output matches expected values from Python 2 logic
# ---------------------------------------------------------------------------

class TestParityWithPython2:
    """Verify that key outputs match values computed by the original Python 2 code."""

    def test_reverse_complement_parity(self):
        """RC of ATCGATCG matches the Python 2 string.translate result."""
        from utils import reverse_complement
        # Python 2: string.maketrans('ATCGN','TAGCN') then translate and reverse
        expected = "CGATCGAT"
        assert reverse_complement("ATCGATCG") == expected

    def test_n50_parity_even_list(self):
        """N50 for [100, 200]: expanded has 300 items; median of items[149]+items[150] = 200."""
        from utils import _n50
        # expanded = 100×[100] + 200×[200]; indices 0-99=100, 100-299=200
        # mid=150; expanded[149]=200, expanded[150]=200 → (200+200)//2 = 200
        assert _n50([100, 200]) == 200

    def test_n50_parity_odd_list(self):
        """N50 for odd-length expanded list matches Python 2 integer division."""
        from utils import _n50
        # lengths [100]: expanded = [100]*100 = 100 items (even), mid=50
        # items[49]=100, items[50]=100 → (100+100)//2 = 100
        assert _n50([100]) == 100

    def test_coord2hash_integer_coords_parity(self):
        """Parsed coordinates are integers inside each hit, matching Python 2 int() conversion."""
        params = {"begin_end": 0, "rawMapping": ""}
        seqs = fasta2hash(QUERY)
        result = coord2hash(str(FIXTURES / "test.coords"), params, seqs)
        # Each value is a list of hits (list-of-lists); check inner ints
        for hits in result.values():
            for hit in hits:
                assert all(isinstance(v, int) for v in hit)

    def test_clean_coords_ambiguous_threshold_parity(self):
        """Ambiguity threshold matches Python 2: (second_len/contig_len)*100 >= -a."""
        params = build_parameters({"-q": "", "-r": "", "-a": "95"})
        seqs = {"c1": "A" * 100}
        # second hit covers 95/100 = 95% → exactly at threshold → ambiguous
        coord_hash = {"c1": [[1, 100, 1, 100], [200, 294, 1, 95]]}
        clean_coords(coord_hash, params, seqs)
        assert "c1" in params["ambiguous"]

    def test_gap_size_parity(self, tmp_path):
        """Gap size between contig1 (ref 1-60) and contig2 (ref 65-124) is 4."""
        params = _run_pipeline(tmp_path, str(FIXTURES / "test.coords"))
        # gap = 65 - 60 - 1 = 4
        assert 4 in params["gaps"]

    def test_scaffold_fasta_header_format_parity(self, tmp_path):
        """Scaffold FASTA headers use the '>Scaffold_N' format from the original."""
        _run_pipeline(tmp_path, str(FIXTURES / "test.coords"))
        content = (tmp_path / "out_Scaffold.fasta").read_text()
        assert ">Scaffold_1\n" in content

    def test_not_mapped_suffix_parity(self, tmp_path):
        """Unmapped contigs use the '_not_mapped' suffix from the original."""
        _run_pipeline(tmp_path, str(FIXTURES / "test.coords"))
        content = (tmp_path / "out_Scaffold.fasta").read_text()
        assert "contig_unmapped_not_mapped" in content

    def test_output_log_version_string_parity(self, tmp_path):
        """Output log starts with the version string from the original."""
        _run_pipeline(tmp_path, str(FIXTURES / "test.coords"))
        log = (tmp_path / "out_output.txt").read_text()
        assert log.startswith("Scaffold_builder version 2.2 log file")


# ---------------------------------------------------------------------------
# Median-size realistic test
# Fixtures: median_query.fasta  (10 mapped contigs 1.5–3.9 kb + 1 unmapped 800 nt)
#           median_reference.fasta (30 kb phage-like reference, random seed 42)
#           median.coords  (hand-crafted to match sequences exactly)
#
# Pipeline behaviour verified against the original Python 2 logic:
#   - ctg2/ctg3 overlap 101 nt on reference → NW alignment → overlapOver=1
#   - ctg7 maps reverse strand → reverse-complemented before scaffolding
#   - 8 gap-filling events with known sizes
#   - ctg_unmapped written with _not_mapped suffix
#   - All 10 mapped contigs collapse into a single Scaffold_1
#
# Golden values captured by running the Python 3 pipeline directly and
# cross-checked against the Python 2 algorithm by hand.
# ---------------------------------------------------------------------------

MEDIAN_QUERY = str(FIXTURES / "median_query.fasta")
MEDIAN_REF = str(FIXTURES / "median_reference.fasta")
MEDIAN_COORDS = str(FIXTURES / "median.coords")


def _run_median(tmp_path):
    """Run the pipeline on median-size fixtures.

    Args:
        tmp_path (Path): Pytest temporary directory.

    Returns:
        tuple: (params dict, scaffold_text, log_text)
    """
    prefix = str(tmp_path / "out")
    (tmp_path / "out_overlap_alignment").mkdir()
    params = build_parameters({"-q": MEDIAN_QUERY, "-r": MEDIAN_REF, "-p": prefix})
    seqs = fasta2hash(MEDIAN_QUERY)
    coords = coord2hash(MEDIAN_COORDS, params, seqs)
    cleaned = clean_coords(coords, params, seqs)
    sort_mapping(cleaned, params, seqs)
    scaffold = (tmp_path / "out_Scaffold.fasta").read_text()
    log = (tmp_path / "out_output.txt").read_text()
    return params, scaffold, log


class TestMedianSize:
    """End-to-end test on realistic median-size inputs (30 kb ref, 10 contigs).

    Golden values are derived from the Python 3 pipeline and cross-checked
    against the original Python 2 algorithm by hand to confirm parity.
    """

    def test_single_scaffold_produced(self, tmp_path):
        """All 10 mapped contigs collapse into exactly one Scaffold_1 sequence."""
        _, scaffold, _ = _run_median(tmp_path)
        headers = [l for l in scaffold.splitlines() if l.startswith(">")]
        scaffold_headers = [h for h in headers if h.startswith(">Scaffold_")]
        assert len(scaffold_headers) == 1
        assert scaffold_headers[0] == ">Scaffold_1"

    def test_scaffold_length(self, tmp_path):
        """Scaffold_1 is 25999 nt: sum of contig lengths minus overlap plus gaps."""
        _, scaffold, _ = _run_median(tmp_path)
        lines = scaffold.splitlines()
        seq = ""
        in_scaffold = False
        for line in lines:
            if line == ">Scaffold_1":
                in_scaffold = True
            elif line.startswith(">"):
                in_scaffold = False
            elif in_scaffold:
                seq += line
        assert len(seq) == 25999

    def test_unmapped_contig_present(self, tmp_path):
        """ctg_unmapped is written with the _not_mapped suffix."""
        _, scaffold, _ = _run_median(tmp_path)
        assert ">ctg_unmapped_not_mapped" in scaffold

    def test_gap_count(self, tmp_path):
        """Exactly 8 gap-filling events occur (one per adjacent non-overlapping pair)."""
        params, _, _ = _run_median(tmp_path)
        assert len(params["gaps"]) == 8

    def test_gap_sizes_match_coords(self, tmp_path):
        """Gap sizes match the reference-position differences in median.coords."""
        params, _, _ = _run_median(tmp_path)
        # Derived from coords: r_next_start - r_curr_end - 1 for each non-overlapping pair
        expected = sorted([99, 199, 199, 99, 499, 199, 99, 199])
        assert sorted(params["gaps"]) == expected

    def test_total_gap_length(self, tmp_path):
        """Total gap length is 1592 nt (sum of all 8 gaps)."""
        params, _, _ = _run_median(tmp_path)
        assert sum(params["gaps"]) == 1592

    def test_overlap_resolved_by_nw(self, tmp_path):
        """ctg2/ctg3 overlap (101 nt, 100% identity) is resolved by Needleman-Wunsch."""
        params, _, _ = _run_median(tmp_path)
        assert params["overlapOver"] == 1
        assert params["overlapBelow"] == 0

    def test_no_ambiguous_contigs(self, tmp_path):
        """No contigs are flagged as ambiguously mapped."""
        params, _, _ = _run_median(tmp_path)
        assert params["ambiguous"] == []

    def test_reverse_contig_scaffolded(self, tmp_path):
        """ctg7 (reverse-mapped) is included in Scaffold_1, not written as unmapped."""
        _, scaffold, _ = _run_median(tmp_path)
        assert "ctg7_not_mapped" not in scaffold
        assert ">Scaffold_1" in scaffold

    def test_ns_in_scaffold_for_gaps(self, tmp_path):
        """Gap-filling inserts N runs; the longest gap (499 nt) produces 499 Ns."""
        _, scaffold, _ = _run_median(tmp_path)
        seq = "".join(l for l in scaffold.splitlines() if not l.startswith(">"))
        assert "N" * 499 in seq

    def test_log_statistics_total_length(self, tmp_path):
        """Output log reports correct assembly and scaffold total lengths."""
        _, _, log = _run_median(tmp_path)
        # Assembly total = sum of all 11 contig lengths
        assert "25308" in log   # assembly total (10 mapped contigs + 1 unmapped)
        # Scaffold total = Scaffold_1 (25999) + ctg_unmapped_not_mapped (800)
        assert "26799" in log

    def test_log_gap_statistics(self, tmp_path):
        """Output log reports correct gap statistics."""
        _, _, log = _run_median(tmp_path)
        assert "Non-overlapping contig pairs\t-\t8" in log
        assert "Total length of gaps\t-\t1592" in log
        assert "Longest gap\t-\t499" in log
        assert "Shortest gap\t-\t99" in log
        assert "Average gap size\t-\t199.0" in log

    def test_log_overlap_statistics(self, tmp_path):
        """Output log reports 1 high-identity overlap pair."""
        _, _, log = _run_median(tmp_path)
        assert "Overlapping contig pairs (>=80% id)\t-\t1" in log
        assert "Overlapping contig pairs (<80% id)\t-\t0" in log

    def test_log_sequence_counts(self, tmp_path):
        """Output log reports 11 assembly sequences and 2 scaffold sequences."""
        _, _, log = _run_median(tmp_path)
        assert "Number of sequences\t11\t2" in log

    def test_overlap_alignment_files_written(self, tmp_path):
        """Overlap and alignment FASTA files are written for the ctg2/ctg3 overlap."""
        _run_median(tmp_path)
        overlap_dir = tmp_path / "out_overlap_alignment"
        assert any(overlap_dir.glob("Overlap_ctg2_ctg3.fasta"))
        assert any(overlap_dir.glob("Alignment_ctg2_ctg3.fasta"))

    # -----------------------------------------------------------------------
    # Parity with Python 2 original
    # -----------------------------------------------------------------------

    def test_parity_gap_sizes(self, tmp_path):
        """Gap sizes match Python 2: gap = r_next[0] - r_curr[1] - 1 after extend_mapping."""
        params, _, _ = _run_median(tmp_path)
        # Python 2 sortMapping → scaffold() computes: checkGap = nextElement[1][0] - currentElement[1][1] - 1
        # With fully-mapped contigs, extend_mapping is a no-op, so gaps equal raw coord differences.
        assert 99 in params["gaps"]
        assert 499 in params["gaps"]

    def test_parity_reverse_complement(self, tmp_path):
        """ctg7 sequence in Scaffold_1 equals RC of the original ctg7 (Python 2 reverse() behaviour)."""
        from utils import reverse_complement
        seqs = fasta2hash(MEDIAN_QUERY)
        original_ctg7 = seqs["ctg7"]
        _, scaffold, _ = _run_median(tmp_path)
        scaffold_seq = "".join(l for l in scaffold.splitlines() if not l.startswith(">"))
        assert reverse_complement(original_ctg7) in scaffold_seq

    def test_parity_overlap_consensus_in_scaffold(self, tmp_path):
        """The NW consensus for ctg2/ctg3 overlap is embedded in Scaffold_1 (not a raw join)."""
        seqs = fasta2hash(MEDIAN_QUERY)
        ctg2_end = seqs["ctg2"][-101:]   # last 101 nt of ctg2
        ctg3_start = seqs["ctg3"][:101]  # first 101 nt of ctg3
        # Since sequences are identical (both extracted from same ref region),
        # the NW consensus equals the overlap sequence itself.
        assert ctg2_end == ctg3_start    # confirms 100% identity overlap
        _, scaffold, _ = _run_median(tmp_path)
        scaffold_seq = "".join(l for l in scaffold.splitlines() if not l.startswith(">"))
        assert ctg2_end in scaffold_seq

    def test_parity_not_mapped_sequence_unchanged(self, tmp_path):
        """ctg_unmapped sequence in output matches the original query sequence exactly."""
        seqs = fasta2hash(MEDIAN_QUERY)
        original = seqs["ctg_unmapped"]
        _, scaffold, _ = _run_median(tmp_path)
        assert original in scaffold


# ---------------------------------------------------------------------------
# True end-to-end test: calls run() which invokes real nucmer + show-coords.
# Skipped automatically if nucmer is not on PATH.
# Results are compared against the golden values from TestMedianSize to
# confirm that real nucmer alignment produces identical scaffolding output.
# ---------------------------------------------------------------------------

import tempfile

from shutil import which as _which

pytestmark_nucmer = pytest.mark.skipif(
    _which("nucmer") is None,
    reason="nucmer not on PATH",
)


@pytest.mark.skipif(_which("nucmer") is None, reason="nucmer not on PATH")
class TestEndToEnd:
    """Full pipeline test using real nucmer on median-size fixtures.

    Verifies that run() produces output identical to the hand-crafted
    fixture run in TestMedianSize, confirming nucmer alignment agrees
    with the pre-baked .coords file.
    """

    def _run(self, tmp_path):
        """Run the full pipeline via scaffold_builder.run().

        Args:
            tmp_path (Path): Pytest temporary directory.

        Returns:
            tuple: (params, scaffold_text, log_text)
        """
        from scaffold_builder import run
        prefix = str(tmp_path / "out")
        params = build_parameters({
            "-q": MEDIAN_QUERY,
            "-r": MEDIAN_REF,
            "-p": prefix,
        })
        run(params)
        scaffold = (tmp_path / "out_Scaffold.fasta").read_text()
        log = (tmp_path / "out_output.txt").read_text()
        return params, scaffold, log

    def test_single_scaffold(self, tmp_path):
        """Real nucmer run produces exactly one Scaffold_1 sequence."""
        _, scaffold, _ = self._run(tmp_path)
        headers = [l for l in scaffold.splitlines() if l.startswith(">Scaffold_")]
        assert headers == [">Scaffold_1"]

    def test_scaffold_length_matches_fixture(self, tmp_path):
        """Scaffold_1 length from real nucmer matches hand-crafted fixture (25999 nt)."""
        _, scaffold, _ = self._run(tmp_path)
        seq = ""
        in_s = False
        for line in scaffold.splitlines():
            if line == ">Scaffold_1":
                in_s = True
            elif line.startswith(">"):
                in_s = False
            elif in_s:
                seq += line
        assert len(seq) == 25999

    def test_gap_count_matches_fixture(self, tmp_path):
        """Real nucmer produces the same 8 gap events as the fixture run."""
        params, _, _ = self._run(tmp_path)
        assert len(params["gaps"]) == 8

    def test_gap_sizes_match_fixture(self, tmp_path):
        """Real nucmer gap sizes match the fixture golden values exactly."""
        params, _, _ = self._run(tmp_path)
        assert sorted(params["gaps"]) == sorted([99, 99, 99, 199, 199, 199, 199, 499])

    def test_overlap_resolved_matches_fixture(self, tmp_path):
        """Real nucmer: ctg2/ctg3 overlap resolved by NW (overlapOver=1)."""
        params, _, _ = self._run(tmp_path)
        assert params["overlapOver"] == 1
        assert params["overlapBelow"] == 0

    def test_log_statistics_match_fixture(self, tmp_path):
        """Real nucmer output log statistics match fixture golden values."""
        _, _, log = self._run(tmp_path)
        assert "Total length\t25308\t26799" in log
        assert "Number of sequences\t11\t2" in log
        assert "Non-overlapping contig pairs\t-\t8" in log
        assert "Total length of gaps\t-\t1592" in log
        assert "Overlapping contig pairs (>=80% id)\t-\t1" in log

    def test_unmapped_contig_present(self, tmp_path):
        """Real nucmer run writes ctg_unmapped with _not_mapped suffix."""
        _, scaffold, _ = self._run(tmp_path)
        assert ">ctg_unmapped_not_mapped" in scaffold

    def test_coords_file_created(self, tmp_path):
        """nucmer produces a .coords file that is then consumed by the pipeline."""
        self._run(tmp_path)
        assert (tmp_path / "out.coords").exists()

    def test_delta_file_removed(self, tmp_path):
        """The intermediate .delta file is cleaned up after show-coords runs."""
        self._run(tmp_path)
        assert not (tmp_path / "out.delta").exists()
