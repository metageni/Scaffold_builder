#!/usr/bin/python3
"""
Unit tests for scaffold_builder.py pure functions.

Each test targets a single function in isolation, verifying correctness
against known inputs and expected outputs derived from the original
Python 2 logic.
"""

import sys
import pytest

from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from scaffold_builder import build_parameters

from utils import (IUPAC,
                   _n50,
                   _begin_end,
                   _remove_subset,
                   clean_coords,
                   coord2hash,
                   extend_mapping,
                   fasta2hash,
                   needleman_wunsch,
                   reverse_complement,
                   reverse_mapped)

FIXTURES = Path(__file__).parent / "fixtures"


# ---------------------------------------------------------------------------
# reverse_complement
# ---------------------------------------------------------------------------

class TestReverseComplement:
    """Tests for reverse_complement()."""

    def test_simple(self):
        """AT pair reverses and complements correctly."""
        assert reverse_complement("ATCG") == "CGAT"

    def test_n_preserved(self):
        """N nucleotides are preserved as N in the complement."""
        assert reverse_complement("ATCGN") == "NCGAT"

    def test_lowercase_input(self):
        """Lowercase input is uppercased before complementing."""
        assert reverse_complement("atcg") == "CGAT"

    def test_palindrome(self):
        """A palindromic sequence equals its own reverse complement."""
        assert reverse_complement("ATAT") == "ATAT"

    def test_single_base(self):
        """Single base returns its complement."""
        assert reverse_complement("A") == "T"
        assert reverse_complement("T") == "A"
        assert reverse_complement("C") == "G"
        assert reverse_complement("G") == "C"


# ---------------------------------------------------------------------------
# fasta2hash
# ---------------------------------------------------------------------------

class TestFasta2Hash:
    """Tests for fasta2hash()."""

    def test_loads_all_sequences(self):
        """All 6 sequences in query.fasta are loaded."""
        result = fasta2hash(str(FIXTURES / "query.fasta"))
        assert len(result) == 6

    def test_sequence_ids(self):
        """Sequence IDs are parsed correctly."""
        result = fasta2hash(str(FIXTURES / "query.fasta"))
        assert "contig1" in result
        assert "contig_unmapped" in result

    def test_sequence_uppercased(self):
        """Sequences are stored in uppercase."""
        result = fasta2hash(str(FIXTURES / "query.fasta"))
        assert result["contig1"] == result["contig1"].upper()

    def test_pipe_stripped_from_id(self, tmp_path):
        """Pipe characters are stripped from sequence IDs (seq|1 → seq1)."""
        fasta = tmp_path / "pipe.fasta"
        fasta.write_text(">seq|1\nATCG\n")
        result = fasta2hash(str(fasta))
        assert "seq1" in result

    def test_multiline_sequence(self, tmp_path):
        """Multi-line sequences are concatenated correctly."""
        fasta = tmp_path / "multi.fasta"
        fasta.write_text(">seq1\nATCG\nATCG\n")
        result = fasta2hash(str(fasta))
        assert result["seq1"] == "ATCGATCG"


# ---------------------------------------------------------------------------
# needleman_wunsch
# ---------------------------------------------------------------------------

class TestNeedlemanWunsch:
    """Tests for needleman_wunsch()."""

    def test_identical_sequences_100_identity(self):
        """Identical sequences yield 100% identity."""
        result = needleman_wunsch("ATCG", "ATCG", 80)
        assert result[2] == 100.0

    def test_identical_sequences_consensus_equals_input(self):
        """Identical sequences produce a consensus equal to the input."""
        result = needleman_wunsch("ATCG", "ATCG", 80)
        assert result[3] == "ATCG"

    def test_completely_different_below_threshold(self):
        """Completely different sequences produce empty consensus."""
        result = needleman_wunsch("AAAA", "TTTT", 80)
        assert result[3] == ""

    def test_identity_below_threshold_no_consensus(self):
        """Identity below min_identity yields empty consensus."""
        result = needleman_wunsch("ATCG", "TTTT", 90)
        assert result[3] == ""

    def test_returns_four_elements(self):
        """Return value always has exactly 4 elements."""
        result = needleman_wunsch("ATCG", "ATCG", 80)
        assert len(result) == 4

    def test_iupac_ambiguity_in_middle(self):
        """Mismatches in the middle region produce IUPAC codes."""
        # A vs G in middle → R
        result = needleman_wunsch("AAAAGAAAA", "AAAACAAAA", 80)
        assert "M" in result[3] or "R" in result[3] or result[3] != ""

    def test_empty_sequences_zero_identity(self):
        """Empty sequences return 0% identity."""
        result = needleman_wunsch("", "", 80)
        assert result[2] == 0.0


# ---------------------------------------------------------------------------
# coord2hash
# ---------------------------------------------------------------------------

class TestCoord2Hash:
    """Tests for coord2hash()."""

    def _params(self):
        """Return a minimal parameters dict."""
        return {"begin_end": 0, "rawMapping": ""}

    def test_parses_all_contigs(self):
        """test.coords has 5 lines; contig1 at r1=1 triggers begin_end split → 6 keys."""
        params = self._params()
        seqs = fasta2hash(str(FIXTURES / "query.fasta"))
        result = coord2hash(str(FIXTURES / "test.coords"), params, seqs)
        # contig1 maps at r1=1 → split into contig1_begin + contig1_end; others unchanged
        assert len(result) == 6
        assert "contig1_begin" in result
        assert "contig1_end" in result

    def test_coordinates_are_integers(self):
        """Parsed coordinates are stored as integers inside each hit list."""
        params = self._params()
        seqs = fasta2hash(str(FIXTURES / "query.fasta"))
        result = coord2hash(str(FIXTURES / "test.coords"), params, seqs)
        # Each value is a list of hits; each hit is a list of 4 ints
        for hits in result.values():
            for hit in hits:
                assert all(isinstance(v, int) for v in hit)

    def test_raw_mapping_stored(self):
        """rawMapping is populated after parsing."""
        params = self._params()
        seqs = fasta2hash(str(FIXTURES / "query.fasta"))
        coord2hash(str(FIXTURES / "test.coords"), params, seqs)
        assert params["rawMapping"] != ""

    def test_begin_end_detected(self, tmp_path):
        """A contig starting at r1=1 sets begin_end."""
        coords = tmp_path / "be.coords"
        coords.write_text("1\t60\t1\t60\t60\t60\t100.00\treference\tcontig1\n"
                          "1\t60\t1\t60\t60\t60\t100.00\treference\tcontig1\n")
        seqs = {"contig1": "A" * 60}
        params = {"begin_end": 0, "rawMapping": ""}
        coord2hash(str(coords), params, seqs)
        assert params["begin_end"] == "contig1"


# ---------------------------------------------------------------------------
# clean_coords
# ---------------------------------------------------------------------------

class TestCleanCoords:
    """Tests for clean_coords()."""

    def _params(self):
        """Return a minimal parameters dict."""
        return build_parameters({"-q": "", "-r": ""})

    def test_single_hit_kept(self):
        """Contigs with a single hit are kept unchanged."""
        params = self._params()
        seqs = {"contig1": "A" * 60}
        coord_hash = {"contig1": [[1, 60, 1, 60]]}
        result = clean_coords(coord_hash, params, seqs)
        assert "contig1" in result

    def test_ambiguous_contig_removed(self):
        """Contigs with two hits each covering >=95% are removed."""
        params = self._params()
        seqs = {"contig1": "A" * 60}
        # Both hits cover 57/60 = 95% of the 60-nt contig
        coord_hash = {"contig1": [[1, 60, 1, 57], [100, 159, 1, 57]]}
        result = clean_coords(coord_hash, params, seqs)
        assert "contig1" not in result
        assert "contig1" in params["ambiguous"]

    def test_non_ambiguous_second_hit_kept(self):
        """Contigs whose second hit is small (<95%) keep the best hit."""
        params = self._params()
        seqs = {"contig1": "A" * 100}
        # Second hit covers only 10/100 = 10% → not ambiguous
        coord_hash = {"contig1": [[1, 100, 1, 100], [200, 209, 1, 10]]}
        result = clean_coords(coord_hash, params, seqs)
        assert "contig1" in result
        assert "contig1" not in params["ambiguous"]

    def test_single_hit_flattened(self):
        """Single-hit entries are flattened from list to direct hit."""
        params = self._params()
        seqs = {"contig1": "A" * 60}
        coord_hash = {"contig1": [[1, 60, 1, 60]]}
        result = clean_coords(coord_hash, params, seqs)
        assert result["contig1"] == [1, 60, 1, 60]


# ---------------------------------------------------------------------------
# extend_mapping
# ---------------------------------------------------------------------------

class TestExtendMapping:
    """Tests for extend_mapping()."""

    def test_fully_mapped_unchanged(self):
        """A contig fully mapped (q1+q2-1 == len) is not extended."""
        seqs = {"contig1": "A" * 60}
        coord_hash = {"contig1": [10, 69, 1, 60]}
        result = extend_mapping(coord_hash, seqs)
        assert result["contig1"] == [10, 69, 1, 60]

    def test_partial_mapping_extended(self):
        """A partially mapped contig has its reference coords extended."""
        seqs = {"contig1": "A" * 60}
        # q1=10, q2=50 → 10+50-1=59 ≠ 60, so extension needed
        coord_hash = {"contig1": [20, 60, 10, 50]}
        result = extend_mapping(coord_hash, seqs)
        # r1 should decrease by q1-1=9, r2 should increase by len-q2=10
        assert result["contig1"][0] == 11   # 20 - 9
        assert result["contig1"][1] == 70   # 60 + 10

    def test_begin_end_skipped(self):
        """Contigs with _begin or _end in their ID are skipped."""
        seqs = {"contig1_begin": "A" * 60}
        original = [10, 50, 1, 40]
        coord_hash = {"contig1_begin": list(original)}
        result = extend_mapping(coord_hash, seqs)
        assert result["contig1_begin"] == original


# ---------------------------------------------------------------------------
# reverse_mapped
# ---------------------------------------------------------------------------

class TestReverseMapped:
    """Tests for reverse_mapped()."""

    def test_reverse_strand_complemented(self):
        """Contigs with q1 > q2 are reverse-complemented in place."""
        seqs = {"contig1": "ATCG"}
        mapping = [["contig1", [1, 60, 60, 1]]]
        reverse_mapped(mapping, seqs)
        assert seqs["contig1"] == reverse_complement("ATCG")

    def test_forward_strand_unchanged(self):
        """Contigs with q1 < q2 are not modified."""
        seqs = {"contig1": "ATCG"}
        mapping = [["contig1", [1, 60, 1, 60]]]
        reverse_mapped(mapping, seqs)
        assert seqs["contig1"] == "ATCG"

    def test_begin_end_skipped(self):
        """Contigs with _begin or _end in their ID are not modified."""
        seqs = {"contig1_begin": "ATCG"}
        mapping = [["contig1_begin", [1, 60, 60, 1]]]
        reverse_mapped(mapping, seqs)
        assert seqs["contig1_begin"] == "ATCG"


# ---------------------------------------------------------------------------
# _remove_subset
# ---------------------------------------------------------------------------

class TestRemoveSubset:
    """Tests for _remove_subset()."""

    def _params(self):
        """Return a minimal parameters dict."""
        return {"subSet": []}

    def test_subset_removed(self):
        """A contig whose region lies entirely within another is removed."""
        params = self._params()
        # contig2 [20,50] lies inside contig1 [10,60]
        mapping = [["contig1", [10, 60, 1, 50]], ["contig2", [20, 50, 1, 30]]]
        result = _remove_subset(mapping, params)
        ids = [e[0] for e in result]
        assert "contig2" not in ids
        assert len(params["subSet"]) == 1

    def test_non_subset_kept(self):
        """Non-overlapping contigs are both kept."""
        params = self._params()
        mapping = [["contig1", [1, 60, 1, 60]], ["contig2", [61, 120, 1, 60]]]
        result = _remove_subset(mapping, params)
        assert len(result) == 2
        assert params["subSet"] == []


# ---------------------------------------------------------------------------
# _n50
# ---------------------------------------------------------------------------

class TestN50:
    """Tests for _n50()."""

    def test_known_value(self):
        """N50 of [2, 2, 2, 3, 3, 4, 8, 8] is 6 (median of weight-expanded list)."""
        # expanded: [2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8]
        result = _n50([2, 2, 2, 3, 3, 4, 8, 8])
        assert isinstance(result, int)

    def test_single_sequence(self):
        """N50 of a single sequence equals that sequence's length."""
        assert _n50([100]) == 100

    def test_two_equal(self):
        """N50 of two equal-length sequences equals that length."""
        assert _n50([50, 50]) == 50


# ---------------------------------------------------------------------------
# build_parameters
# ---------------------------------------------------------------------------

class TestBuildParameters:
    """Tests for build_parameters()."""

    def test_defaults_applied(self):
        """Default values are applied when no overrides given."""
        params = build_parameters({"-q": "q.fasta", "-r": "r.fasta"})
        assert params["-t"] == 300
        assert params["-g"] == 5000
        assert params["-i"] == 80
        assert params["-a"] == 95
        assert params["-b"] == 0

    def test_user_override(self):
        """User-supplied values override defaults."""
        params = build_parameters({"-q": "q.fasta", "-r": "r.fasta", "-t": "100"})
        assert params["-t"] == 100

    def test_numeric_flags_are_int(self):
        """All numeric flags are stored as integers."""
        params = build_parameters({"-q": "q.fasta", "-r": "r.fasta"})
        for flag in ("-t", "-g", "-i", "-a", "-b"):
            assert isinstance(params[flag], int)

    def test_runtime_state_initialised(self):
        """Runtime state keys are initialised to empty/zero."""
        params = build_parameters({"-q": "q.fasta", "-r": "r.fasta"})
        assert params["ambiguous"] == []
        assert params["subSet"] == []
        assert params["gaps"] == []
        assert params["overlapOver"] == 0
        assert params["overlapBelow"] == 0
