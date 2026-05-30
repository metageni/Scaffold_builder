# Scaffold_builder Project

## Purpose
Python 3 rewrite of Scaffold_builder v2.2 — a bioinformatics tool that scaffolds
de novo assembled contigs against a reference genome using MUMmer (nucmer) alignments.

## Key Algorithms
- **Needleman-Wunsch** alignment for overlap resolution between adjacent contigs
- **IUPAC ambiguity codes** for consensus generation at mismatched positions
- **Ambiguity detection**: contigs mapping ≥ `-a`% of their length to 2+ locations are excluded
- **Gap filling**: non-overlapping adjacent contigs joined with N-padding
- **Circular reference**: contigs mapping at r1=1 are split into `_begin`/`_end` entries

## Parameters (DEFAULTS in utils.py)
| Flag | Default | Meaning |
|------|---------|---------|
| `-a` | 95 | Ambiguous mapping threshold (%) |
| `-i` | 80 | Minimum identity for merging overlaps (%) |
| `-t` | 300 | Terminus length to align (nt) |
| `-g` | 5000 | Maximum gap length before scaffold break (nt) |
| `-b` | 0 | Rearrangement behaviour: 0=end-to-end, 1=new scaffold |

## Pipeline Order
1. `run_nucmer()` — runs nucmer + show-coords, writes `<prefix>.coords`
2. `fasta2hash()` — loads query contigs into memory
3. `coord2hash()` — parses .coords; detects circular begin/end
4. `clean_coords()` — removes ambiguously mapped contigs
5. `sort_mapping()` → `extend_mapping` → `reverse_mapped` → `_remove_subset` → `scaffold_sequences`
6. `write_output()` — writes log + statistics

## Output Files
- `<prefix>_Scaffold.fasta` — scaffolded sequences + unmapped/ambiguous contigs
- `<prefix>_output.txt` — log with statistics
- `<prefix>_overlap_alignment/` — per-pair Overlap_*.fasta and Alignment_*.fasta

## Testing
- Run: `python3 -m pytest tests/ -v`
- 112 tests total: 40 unit (test_unit.py) + 72 integration (test_integration.py)
- Fixtures in `tests/fixtures/`: query.fasta, reference.fasta, test.coords, ambiguous.coords, reverse.coords, median_query.fasta, median_reference.fasta, median.coords, large_query.fasta, large_reference.fasta, large.coords, large2_query.fasta, large2_reference.fasta, large3_query.fasta, large3_reference.fasta
- Integration tests bypass nucmer entirely using pre-baked .coords fixtures
- nucmer IS available in the dev environment; end-to-end tests (TestLargeEndToEnd, TestLarge2EndToEnd, TestLarge3EndToEnd) call `run_nucmer()` directly

## Packaging
- `pyproject.toml` defines the package; entry point: `scaffold_builder = "cli:main"`
- Runtime dependency: `click>=8.0`
- Build: `pip install .`

## Parity Requirement
All logic must produce output identical to the original Python 2 scaffold_builder.py.
The `TestParityWithPython2` class in test_integration.py documents known equivalences.
