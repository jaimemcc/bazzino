from __future__ import annotations

import argparse
import sys
from pathlib import Path

import jupytext

from notebook_workflow_config import collect_tracked_notebooks
from notebook_workflow_config import file_digest
from notebook_workflow_config import load_workflow_state
from notebook_workflow_config import load_managed_roots
from notebook_workflow_config import record_workflow_state
from notebook_workflow_config import save_workflow_state
from notebook_workflow_config import state_key
from notebook_workflow_config import resolve_tracked_notebook_arg
from notebook_workflow_config import source_path_for_tracked


ROOT = Path(__file__).resolve().parents[2]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Regenerate source .ipynb notebooks from tracked .py files.")
    parser.add_argument(
        "notebook",
        nargs="?",
        help="Specific notebook to regenerate (e.g., 'starter_notebook.py'). If omitted, regenerates all notebooks.",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite source files when both sides changed.")
    args = parser.parse_args(argv)

    try:
        managed_roots = load_managed_roots(ROOT)
    except ValueError as exc:
        print(f"Notebook workflow config error: {exc}", file=sys.stderr)
        return 2

    if args.notebook:
        tracked_notebook = resolve_tracked_notebook_arg(args.notebook, managed_roots, ROOT)
        if tracked_notebook is None:
            print(
                f"Notebook not found or ambiguous: {args.notebook}. "
                "Try a path relative to repository root, such as 'feature1/notebooks/text/example.py'.",
                file=sys.stderr,
            )
            return 1
        notebook_entries = [
            (managed_root, tracked_notebook)
            for managed_root, candidate in collect_tracked_notebooks(managed_roots)
            if candidate == tracked_notebook
        ]
    else:
        notebook_entries = collect_tracked_notebooks(managed_roots)
        if not notebook_entries:
            print("No tracked notebooks found under configured managed roots.")
            return 0

    try:
        state = load_workflow_state(ROOT)
    except ValueError as exc:
        print(f"Notebook workflow state error: {exc}", file=sys.stderr)
        return 2

    changed = False
    for managed_root, tracked_notebook in notebook_entries:
        source_notebook = source_path_for_tracked(tracked_notebook, managed_root)
        source_notebook.parent.mkdir(parents=True, exist_ok=True)

        notebook_object = jupytext.read(tracked_notebook, fmt="py:percent")
        regenerated_text = jupytext.writes(notebook_object, fmt="ipynb")
        source_key = state_key(source_notebook, ROOT)
        previous = state.get("notebooks", {}).get(source_key)
        tracked_changed = not previous or previous.get("tracked_sha256") != file_digest(tracked_notebook)
        source_changed = source_notebook.exists() and (
            not previous or previous.get("source_sha256") != file_digest(source_notebook)
        )
        content_changed = source_notebook.exists() and source_notebook.read_text(encoding="utf-8") != regenerated_text

        if content_changed and source_changed and tracked_changed and not args.force:
            print(
                f"Regeneration conflict for {source_notebook.relative_to(ROOT)}: both the .ipynb and .py changed "
                f"since the last recorded sync. Resolve manually, or run 'pixi run regenerate-notebooks --force' "
                f"to make the .py authoritative.",
                file=sys.stderr,
            )
            return 1
        if content_changed and source_changed and not args.force:
            print(
                f"Regeneration refused for {source_notebook.relative_to(ROOT)}: the existing .ipynb has "
                f"changed since the last recorded sync. Use 'pixi run regenerate-notebooks --force' to replace it.",
                file=sys.stderr,
            )
            return 1

        if not source_notebook.exists() or content_changed:
            source_notebook.write_text(regenerated_text, encoding="utf-8")
        record_workflow_state(
            state,
            source_notebook=source_notebook,
            tracked_notebook=tracked_notebook,
            root=ROOT,
            operation="regen",
        )
        changed = True

    if changed:
        save_workflow_state(ROOT, state)

    print("Notebook regeneration complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
