"""Run the test suite and record what it actually reported.

A paper that claims a tool is tested should state a number a reader can reproduce, and that number should
come from running the suite rather than from counting files. This runs pytest, parses its own summary, and
writes the counts alongside every other benchmark result, so the manuscript binds them like any other
value and a claim about coverage cannot drift from the suite.

    python benchmarks/testsuite.py

Output: benchmarks/raw/testsuite.json
"""
from __future__ import annotations
import json, os, subprocess, sys, tempfile, time
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "benchmarks", "raw", "testsuite.json")


def main():
    t0 = time.time()
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    # Counts come from the JUnit report rather than from parsing the terminal output. Under this
    # project's pytest configuration a fully clean `-q` run prints only the progress dots and no summary
    # line at all, so a parser reports zero tests for a green suite - which is the one failure a number
    # in a manuscript must not have. The XML carries the counts as attributes.
    with tempfile.TemporaryDirectory() as tmp:
        xml = os.path.join(tmp, "junit.xml")
        p = subprocess.run([sys.executable, "-m", "pytest", "tests", "-q", "--no-header",
                            f"--junitxml={xml}"],
                           cwd=ROOT, capture_output=True, text=True, encoding="utf-8",
                           errors="replace", env=env)
        if not os.path.exists(xml):
            print("pytest produced no JUnit report; not writing a result")
            print((p.stdout or "")[-2000:])
            return 1
        root = ET.parse(xml).getroot()
        suite = root if root.tag == "testsuite" else root.find("testsuite")
        a = suite.attrib
        total = int(a.get("tests", 0))
        failures = int(a.get("failures", 0))
        errors = int(a.get("errors", 0))
        skipped = int(a.get("skipped", 0))

    files = sorted(f for f in os.listdir(os.path.join(ROOT, "tests")) if f.startswith("test_")
                   and f.endswith(".py"))
    result = {
        "exit_code": p.returncode,
        "green": p.returncode == 0 and failures == 0 and errors == 0,
        "counts": {"tests": total, "failures": failures, "errors": errors, "skipped": skipped},
        "passed": total - failures - errors - skipped,
        "failed": failures + errors,
        "skipped": skipped,
        "collected": total,
        "test_files": len(files),
        "elapsed_s": round(time.time() - t0, 1),
        "python": sys.version.split()[0],
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(result, open(OUT, "w", encoding="utf-8"), indent=1, sort_keys=True)
    print(f"exit {p.returncode}")
    print(f"{result['passed']} passed, {result['failed']} failed, {result['skipped']} skipped "
          f"of {result['collected']} collected across {result['test_files']} test files "
          f"in {result['elapsed_s']} s")
    print(f"-> {OUT}")
    return 0 if result["green"] else 1


if __name__ == "__main__":
    sys.exit(main())
