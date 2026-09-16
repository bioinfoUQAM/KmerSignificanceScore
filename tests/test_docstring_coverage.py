"""
Check that every definition of the package carries a docstring.

The contract a reader audits is the one written next to the code. Coverage that holds today
drifts the moment a function is added, and the gap is invisible: nothing fails, the module
imports, and the next reader meets a definition with no statement of what it returns.

This walks the syntax tree rather than searching the text, so a string that merely sits at the
top of a body does not pass for a docstring.

Run directly:
    python tests/test_docstring_coverage.py
or under pytest:
    pytest tests/test_docstring_coverage.py
"""

from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFINITION = (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)


def source_modules() -> list[Path]:
    """The modules the package is audited on: src/*.py and the entry point."""
    return sorted((ROOT / "src").glob("*.py")) + [ROOT / "main.py"]


def count_documented_definitions() -> int:
    """Count the definitions of src/ and main.py, failing if any lacks a docstring."""
    undocumented, definitions = [], 0
    for module in source_modules():
        tree = ast.parse(module.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, DEFINITION):
                definitions += 1
                if not ast.get_docstring(node):
                    undocumented.append(f"{module.name}:{node.lineno}:{node.name}")
    assert definitions, "no definition found; the source tree was not read"
    assert not undocumented, (
        f"{len(undocumented)} definition(s) without a docstring: {', '.join(undocumented)}"
    )
    return definitions


def test_every_definition_has_a_docstring() -> None:
    """pytest wrapper: the counting helper returns a value, a test must not."""
    count_documented_definitions()


if __name__ == "__main__":
    total = count_documented_definitions()
    print(f"every definition carries a docstring    {total} definitions")
