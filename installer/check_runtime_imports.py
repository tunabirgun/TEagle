"""Fail when shipped Python source imports an undeclared runtime package."""
from __future__ import annotations

import argparse
import ast
import re
import sys
from pathlib import Path


IMPORT_BY_DISTRIBUTION = {
    "primer3-py": "primer3",
    "pyhmmer": "pyhmmer",
    "psutil": "psutil",
    "PySide6": "PySide6",
    "openpyxl": "openpyxl",
    "et_xmlfile": "et_xmlfile",
    "certifi": "certifi",
    "pypdf": "pypdf",
}
OPTIONAL_IMPORTS = {"PIL", "RNA"}
TRANSITIVE_PARENT = {"psutil": "pyhmmer", "et_xmlfile": "openpyxl"}


def _requirements(path: Path) -> set[str]:
    declared: set[str] = set()
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        match = re.fullmatch(r"([^=<>!~\s]+)==([^\s]+)", line)
        if not match:
            raise ValueError(f"runtime requirement is not exactly pinned: {line}")
        distribution = match.group(1)
        if distribution in declared:
            raise ValueError(f"duplicate runtime requirement: {distribution}")
        if distribution not in IMPORT_BY_DISTRIBUTION:
            raise ValueError(f"runtime requirement lacks an import mapping: {distribution}")
        declared.add(distribution)
    return declared


def _runtime_files(root: Path) -> list[Path]:
    backend = root / "app" / "backend"
    files = sorted(backend.glob("*.py")) + sorted((backend / "teagle_core").rglob("*.py"))
    files += sorted((root / "app" / "native").glob("*.py"))
    files.append(root / "app" / "teagle.py")
    files.append(root / "installer" / "teagle_native.py")
    return [path for path in files if path.is_file()]


def undeclared_imports(root: Path) -> dict[str, list[str]]:
    declared = _requirements(root / "app" / "backend" / "requirements.txt")
    allowed_external = {IMPORT_BY_DISTRIBUTION[name] for name in declared} | OPTIONAL_IMPORTS
    files = _runtime_files(root)
    internal = {path.stem for path in files} | {"teagle_core"}
    findings: dict[str, list[str]] = {}
    imported_all: set[str] = set()
    for path in files:
        tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(alias.name.split(".", 1)[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                imported.add(node.module.split(".", 1)[0])
        for module in imported:
            imported_all.add(module)
            if module in sys.stdlib_module_names or module in internal or module in allowed_external:
                continue
            findings.setdefault(module, []).append(str(path.relative_to(root)))
    used_distributions = {
        distribution for distribution in declared
        if IMPORT_BY_DISTRIBUTION[distribution] in imported_all
    }
    used_distributions.update(
        distribution for distribution, parent in TRANSITIVE_PARENT.items()
        if distribution in declared and parent in declared
        and IMPORT_BY_DISTRIBUTION[parent] in imported_all
    )
    for distribution, parent in sorted(TRANSITIVE_PARENT.items()):
        if parent in declared and distribution not in declared:
            findings[f"missing transitive requirement {distribution} (required by {parent})"] = [
                "app/backend/requirements.txt"
            ]
    for distribution in sorted(declared - used_distributions):
        findings[f"unused requirement {distribution}"] = ["app/backend/requirements.txt"]
    return findings


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    try:
        findings = undeclared_imports(args.root.resolve())
    except (OSError, SyntaxError, ValueError) as exc:
        print(f"Runtime import parity check failed: {exc}", file=sys.stderr)
        return 2
    if findings:
        for finding, paths in sorted(findings.items()):
            label = (finding if finding.startswith(("missing ", "unused "))
                     else f"undeclared runtime import {finding}")
            print(f"{label}: {', '.join(paths)}", file=sys.stderr)
        return 1
    print("Runtime import parity: exact")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
