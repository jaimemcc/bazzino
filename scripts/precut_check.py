"""Pre-release checks for the publication repo. Run before cutting the new repo.

Usage:
    pixi run python scripts/precut_check.py

Checks tracked files (plus data/assembled/, which ships but may be untracked) for:
  - personal paths / usernames in text content
  - files that must not ship (caches, pickles, .afdesign, panels, notebook text copies)
  - saved notebook errors and autoreload leftovers
  - large files
  - required files and folders
Exits with status 1 if any problem is found.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

FORBIDDEN_PATTERNS = [
    r"jmc010",
    r"D:[/\\]+TestData",
    r"C:[/\\]+Users[/\\]+",
]
FORBIDDEN_FILES = [
    r"\.pickle$", r"\.pkl$", r"\.afdesign$", r"_cache_", r"^paper/figs/panels/",
    r"^notebooks/text/", r"^tooling/", r"\.egg-info/", r"\.pyc$", r"^\.pixi/",
]
REQUIRED = [
    "README.md", "LICENSE", "CITATION.cff", "pixi.toml", "pixi.lock",
    "src/README.md", "src/data_io.py", "scripts/download_data.py",
    "data/assembled/manifest.json", "data/assembled/arrays.npz",
]
# plan/ documents the process and legitimately mentions old paths; it is not shipped.
SKIP_CONTENT_PREFIXES = ("plan/", "scripts/precut_check.py", "pixi.lock")
TEXT_SUFFIXES = {".py", ".md", ".ipynb", ".toml", ".cff", ".json", ".csv", ".txt", ".yaml", ".yml", ".cfg"}
LARGE_MB = 25


def tracked_files():
    out = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True).stdout
    files = {Path(p) for p in out.splitlines()}
    assembled = ROOT / "data" / "assembled"
    if assembled.exists():
        files |= {p.relative_to(ROOT) for p in assembled.rglob("*") if p.is_file()}
    return sorted(f for f in files if (ROOT / f).exists())


def main():
    problems = []
    files = tracked_files()

    for rel in files:
        posix = rel.as_posix()
        for pat in FORBIDDEN_FILES:
            if re.search(pat, posix):
                problems.append(f"forbidden file: {posix}")
        if posix.startswith(SKIP_CONTENT_PREFIXES) or rel.suffix.lower() not in TEXT_SUFFIXES:
            continue
        text = (ROOT / rel).read_text(encoding="utf8", errors="ignore")
        for pat in FORBIDDEN_PATTERNS:
            if re.search(pat, text, flags=re.I):
                problems.append(f"personal path /{pat}/ in {posix}")
        if rel.suffix == ".ipynb":
            nb = json.loads(text)
            for i, cell in enumerate(nb["cells"]):
                if cell["cell_type"] != "code":
                    continue
                if any(o.get("output_type") == "error" for o in cell.get("outputs", [])):
                    problems.append(f"saved error output: {posix} cell {i}")
                if any("autoreload" in line for line in cell["source"]):
                    problems.append(f"autoreload left in {posix} cell {i}")
            if not any(c["cell_type"] == "code" and c.get("outputs") for c in nb["cells"]):
                problems.append(f"notebook has no saved outputs: {posix}")

    for req in REQUIRED:
        if not (ROOT / req).exists():
            problems.append(f"missing required: {req}")

    sizes = sorted(((ROOT / f).stat().st_size, f.as_posix()) for f in files)
    total = sum(s for s, _ in sizes)
    print(f"{len(files)} files, {total / 1e6:.1f} MB total. Largest:")
    for size, name in sizes[-8:][::-1]:
        flag = "  <-- large" if size > LARGE_MB * 1e6 else ""
        print(f"  {size / 1e6:7.2f} MB  {name}{flag}")
        if flag:
            problems.append(f"large file (> {LARGE_MB} MB): {name}")

    print()
    if problems:
        print(f"{len(problems)} problem(s):")
        for p in sorted(set(problems)):
            print("  -", p)
        sys.exit(1)
    print("All checks passed.")


if __name__ == "__main__":
    main()
