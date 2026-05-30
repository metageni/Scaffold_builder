#!/usr/bin/env python3
"""
Parity tests: compare Python 3 rewrite output against the original Python 2
scaffold_builder.py, executed via a minimal Python 3 compatibility shim.

The shim patches exactly four Python 2-only constructs:
  - file() → open()
  - string.maketrans → str.maketrans
  - integer division in N50 (len/2 → len//2)
  - print statements → print()

The scaffold FASTA output must be byte-for-byte identical.
The log output must match on all lines except the version string header
and the N50 row (the old code has a float division bug: 2501.0 vs 2501).
"""

import os
import sys
import shutil
import string
import tempfile

import pytest

from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from scaffold_builder import build_parameters

from utils import (clean_coords,
                   coord2hash,
                   fasta2hash,
                   sort_mapping)

FIXTURES = Path(__file__).parent / "fixtures"
OLD_SRC  = Path(__file__).parent.parent.parent / "old" / "Scaffold_builder" / "scaffold_builder.py"
QUERY    = str(FIXTURES / "median_query.fasta")
REF      = str(FIXTURES / "median_reference.fasta")
COORDS   = str(FIXTURES / "median.coords")


def _patch_old_src(src):
    """Apply minimal Python 3 compatibility patches to the old source.

    Args:
        src (str): Raw Python 2 source text.

    Returns:
        str: Patched source runnable under Python 3.
    """
    src = src.replace("file(fasta)", "open(fasta)")
    src = src.replace(
        "complement = string.maketrans('ATCGN', 'TAGCN')",
        "complement = str.maketrans('ATCGN', 'TAGCN')",
    )
    src = src.replace(
        "else:median=finalListLength[len(finalListLength)/2]",
        "else:median=finalListLength[len(finalListLength)//2]",
    )
    src = src.replace('print userParameters\n', 'print(userParameters)\n')
    src = src.replace('print "#Please inform the query file -q"',
                      'print("#Please inform the query file -q")')
    src = src.replace('print "#Please inform the reference file -r"',
                      'print("#Please inform the reference file -r")')
    src = src.replace('print "Scaffolded :)"', 'print("Scaffolded :)")')
    src = src.replace('print helpMsg', 'print(helpMsg)')
    src = src.replace('![](logo/scaffold_builder_logo.png "Logo")', '')
    return src


def _run_old(prefix):
    """Execute the patched old pipeline on the median fixtures.

    Args:
        prefix (str): Output file prefix (directory must exist).

    Returns:
        dict: The parameters dict after the run.
    """
    os.makedirs(prefix + "_overlap_alignment", exist_ok=True)
    shutil.copy(COORDS, prefix + ".coords")

    src = _patch_old_src(OLD_SRC.read_text())
    ns = {"__name__": "__not_main__", "string": string, "os": os, "sys": sys}
    exec(compile(src, "scaffold_builder_old.py", "exec"), ns)

    p = ns["parameters"]
    p["-q"] = QUERY
    p["-r"] = REF
    p["-p"] = prefix
    for flag, val in [("-t", 300), ("-g", 5000), ("-i", 80), ("-a", 95), ("-b", 0)]:
        p[flag] = val

    ns["hashSequences"] = ns["fasta2hash"](QUERY)
    coords = ns["coord2hash"](prefix + ".coords")
    cleaned = ns["cleanCoords"](coords)
    ns["sortMapping"](cleaned)
    return p


def _run_new(prefix):
    """Execute the Python 3 pipeline on the median fixtures.

    Args:
        prefix (str): Output file prefix (directory must exist).

    Returns:
        dict: The parameters dict after the run.
    """
    os.makedirs(prefix + "_overlap_alignment", exist_ok=True)
    shutil.copy(COORDS, prefix + ".coords")

    params = build_parameters({"-q": QUERY, "-r": REF, "-p": prefix})
    seqs = fasta2hash(QUERY)
    coords = coord2hash(prefix + ".coords", params, seqs)
    cleaned = clean_coords(coords, params, seqs)
    sort_mapping(cleaned, params, seqs)
    return params


@pytest.fixture(scope="module")
def both_runs():
    """Run both pipelines once and return their outputs.

    Yields:
        tuple: (old_scaffold, new_scaffold, old_log, new_log, old_params, new_params)
    """
    tmp = tempfile.mkdtemp()
    try:
        old_p = _run_old(os.path.join(tmp, "old"))
        new_p = _run_new(os.path.join(tmp, "new"))
        old_scaffold = open(os.path.join(tmp, "old_Scaffold.fasta")).read()
        new_scaffold = open(os.path.join(tmp, "new_Scaffold.fasta")).read()
        old_log = open(os.path.join(tmp, "old_output.txt")).read()
        new_log = open(os.path.join(tmp, "new_output.txt")).read()
        yield old_scaffold, new_scaffold, old_log, new_log, old_p, new_p
    finally:
        shutil.rmtree(tmp)


@pytest.mark.skipif(not OLD_SRC.exists(), reason="old scaffold_builder.py not found")
class TestParityOldVsNew:
    """Compare Python 3 rewrite output against the patched Python 2 original."""

    def test_scaffold_fasta_identical(self, both_runs):
        """Scaffold FASTA output is byte-for-byte identical between old and new."""
        old_scaffold, new_scaffold, *_ = both_runs
        assert old_scaffold == new_scaffold

    def test_gaps_identical(self, both_runs):
        """Gap list is identical between old and new."""
        *_, old_p, new_p = both_runs
        assert sorted(old_p["gaps"]) == sorted(new_p["gaps"])

    def test_overlap_counts_identical(self, both_runs):
        """overlapOver and overlapBelow counts are identical."""
        *_, old_p, new_p = both_runs
        assert old_p["overlapOver"] == new_p["overlapOver"]
        assert old_p["overlapBelow"] == new_p["overlapBelow"]

    def test_ambiguous_identical(self, both_runs):
        """Ambiguous contig list is identical."""
        *_, old_p, new_p = both_runs
        assert old_p["ambiguous"] == new_p["ambiguous"]

    def test_log_matches_except_version_and_n50(self, both_runs):
        """Log output matches on all lines except version header and N50 row.

        The old code has a Python 2 float-division bug in N50 (returns 2501.0
        instead of 2501). All other log lines must be identical.
        """
        _, _, old_log, new_log, *_ = both_runs
        skip = {"Scaffold_builder version", "N50\t"}
        old_lines = [l for l in old_log.splitlines() if not any(l.startswith(s) for s in skip)]
        new_lines = [l for l in new_log.splitlines() if not any(l.startswith(s) for s in skip)]
        assert old_lines == new_lines
