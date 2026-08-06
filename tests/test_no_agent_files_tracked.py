"""No AI-assistant instruction or configuration file is tracked in this repository.

The project's position is that tooling instructions are not part of the software. They describe how a
particular assistant should behave on this machine, they date faster than anything they describe, and
they have no bearing on what TEagle computes — so they belong beside the checkout, not inside it.

Stating that in `.gitignore` does not enforce it here. `.gitignore` ignores itself (line 23) and is
untracked, so its rules never reach a clone: a contributor working from a fresh clone has no ignore file
at all, and `git add -A` picks up whatever an assistant has written into the working directory. The rule
has to live somewhere that travels with the code.

This is the repository's own idiom applied to itself — prefer a claim the suite can check over one a
reader must verify. It runs `git ls-files`, so it tests the index rather than the disk: the files may sit
in the working directory, where they are useful, and the test only objects when one is staged or
committed.
"""
from __future__ import annotations
import os
import re
import subprocess

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Matched against the basename, case-insensitively. Deliberately broader than the two files that
# prompted it: the next assistant will not necessarily use the same filename, and a guard that only
# knows today's names stops working the moment the tooling changes.
FORBIDDEN = re.compile(
    r"""^(
        agents?\.md            |   # AGENTS.md
        claude\.md             |   # CLAUDE.md
        \.claude.*             |   # .claude/, .claude.json
        \.cursor.*             |   # .cursorrules, .cursor/
        \.aider.*              |   # .aider.conf.yml, .aider.chat.history.md
        \.windsurfrules        |
        \.continuerules        |
        copilot-instructions\.md   |
        \.codeium.*            |
        gemini\.md             |
        \.goosehints
    )$""",
    re.IGNORECASE | re.VERBOSE,
)


def tracked_files():
    r = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        pytest.skip("not a git checkout, or git is unavailable")
    return [line for line in r.stdout.splitlines() if line.strip()]


def test_git_is_actually_reporting_files():
    """Guard the guard. If `git ls-files` returned nothing — not a checkout, git missing, a cwd slip —
    the check below would pass on an empty list and silently stop protecting anything."""
    files = tracked_files()
    assert len(files) > 50, (
        f"git ls-files reported only {len(files)} tracked files; this test cannot conclude anything "
        f"about the index from that, so treat it as broken rather than as passing")


def test_no_ai_assistant_files_are_tracked():
    offenders = sorted(p for p in tracked_files() if FORBIDDEN.match(os.path.basename(p)))
    assert not offenders, (
        "AI-assistant instruction or configuration files are tracked in this repository: "
        + ", ".join(offenders)
        + ". They are not part of the software. Untrack with `git rm --cached <path>`, which leaves the "
          "file in your working directory where it is still read.")
