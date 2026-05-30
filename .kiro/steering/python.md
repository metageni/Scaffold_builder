# Python Coding Standards

## File Structure
- First line: `#!/usr/bin/python3`
- Second: module-level docstring describing the file's purpose
- Last line: blank (POSIX compliance)

## Imports
Sort within each group by **line length (shortest first)**. Groups in order:
1. Native `import` statements
2. *(blank line)*
3. Third-party `import` statements
4. *(blank line)*
5. Native `from ... import` statements
6. *(blank line)*
7. Third-party `from ... import` statements
8. *(blank line)*
9. Local `import` statements
10. *(blank line)*
11. Local `from ... import` statements

Multi-symbol `from` imports use parentheses with each symbol on its own line:
```python
from utils import (clean_coords,
                   coord2hash)
```

## Docstrings
All functions, methods, and classes must have PEP 257 docstrings with Args and Returns sections.

## Logging
Use `logging` instead of `print()`. Reserve `print()` only for explicit CLI stdout output.
Configure `logging.basicConfig()` at the entry point only (`cli.py` / `main.py`).

## CLI Separation
CLI code lives exclusively in `cli.py`. Core logic modules must never contain argparse or sys.argv.
CLI is implemented with `click` (not argparse).

## Module Layout
```
utils.py            — pure/stateless functions only
scaffold_builder.py — pipeline orchestration (build_parameters, run_nucmer, run)
cli.py              — click CLI + logging config
main.py             — entry point: calls cli.main()
tests/
  test_unit.py        — unit tests; import pure functions from utils
  test_integration.py — integration tests; import from scaffold_builder + utils
```
