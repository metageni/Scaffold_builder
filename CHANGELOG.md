# Changelog

## [3.0.0] — 2026-05-30

### Rewrite
Complete rewrite from Python 2 to Python 3. The codebase was restructured into
four modules:

- `utils.py` — pure/stateless functions (alignment, parsing, scaffolding logic)
- `scaffold_builder.py` — pipeline orchestration (`build_parameters`, `run_nucmer`, `run`)
- `cli.py` — command-line interface (click)
- `main.py` — entry point

The rewrite produces **byte-for-byte identical scaffold FASTA output** to v2.2
on all tested inputs (verified by end-to-end comparison against the original
Python 2 binary running in Docker).

### Bug Fixes

- **N50 calculation** — the original code used integer division (`len(list)/2`)
  to index into the sorted length list, which selected the wrong contig instead
  of walking the cumulative-sum threshold. The correct N50 is now computed by
  accumulating lengths in descending order until 50% of total assembly length
  is reached. Example on the large dataset: old reported 176405, correct value
  is 279979.

- **Average length** — the original code performed Python 2 floor division
  (`total/count`), truncating the decimal. Now computed with true float
  division. Example: old reported 12101.0, correct value is 12101.54.

### Other Changes
- Replaced `os.system()` nucmer call with `subprocess.run()` (safer, raises on
  failure instead of silently continuing).
- Replaced `print` statements with `logging` module; log level configurable at
  the CLI entry point.
- Packaged via `pyproject.toml`; installable with `pip install .`.
- Entry point: `scaffold_builder` command (previously required calling
  `python scaffold_builder.py` directly).
- Added test suite: 114 tests covering unit, integration, and end-to-end parity
  against the original Python 2 binary.

---

## [2.2] — 2013

Original Python 2 release by Genivaldo G. Z. Silva.

> Silva GG, Dutilh BE, Matthews TD, Elkins K, Schmieder R, Dinsdale EA,
> Edwards RA. Combining de novo and reference-guided assembly with
> Scaffold_builder. Source Code for Biology and Medicine, 2013.
