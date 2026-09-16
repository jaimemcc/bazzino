from __future__ import annotations

from pathlib import Path

import jupytext

from notebook_workflow_config import collect_source_notebooks
from notebook_workflow_config import file_digest
from notebook_workflow_config import load_workflow_state
from notebook_workflow_config import load_managed_roots
from notebook_workflow_config import record_workflow_state
from notebook_workflow_config import save_workflow_state
from notebook_workflow_config import state_key
from notebook_workflow_config import tracked_path_for_source


ROOT = Path(__file__).resolve().parents[2]


def read_source_notebook(source_notebook: Path):
    return jupytext.reads(source_notebook.read_text(encoding="utf-8-sig"), fmt="ipynb")


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Sync source .ipynb notebooks to tracked .py files.")
    parser.add_argument("--force", action="store_true", help="Overwrite tracked files when both sides changed.")
    args = parser.parse_args(argv)

    try:
        managed_roots = load_managed_roots(ROOT)
    except ValueError as exc:
        print(f"Notebook workflow config error: {exc}")
        return 2

    notebooks = collect_source_notebooks(managed_roots)
    if not notebooks:
        print("No source notebooks found under configured managed roots.")
        return 0

    try:
        state = load_workflow_state(ROOT)
    except ValueError as exc:
        print(f"Notebook workflow state error: {exc}")
        return 2

    pending: list[tuple[Path, Path, str, bool]] = []
    for managed_root, source_notebook in notebooks:
        target_notebook = tracked_path_for_source(source_notebook, managed_root)
        notebook_object = read_source_notebook(source_notebook)
        regenerated_text = jupytext.writes(notebook_object, fmt="py:percent")
        source_key = state_key(source_notebook, ROOT)
        previous = state.get("notebooks", {}).get(source_key)
        source_changed = not previous or previous.get("source_sha256") != file_digest(source_notebook)
        tracked_changed = target_notebook.exists() and (
            not previous or previous.get("tracked_sha256") != file_digest(target_notebook)
        )
        content_changed = target_notebook.exists() and target_notebook.read_text(encoding="utf-8") != regenerated_text

        if content_changed and tracked_changed and source_changed and not args.force:
            print(
                f"Sync conflict for {source_notebook.relative_to(ROOT)}: both the .ipynb and .py changed "
                f"since the last recorded sync. Resolve manually, or run 'pixi run sync-notebooks --force' "
                f"to make the .ipynb authoritative."
            )
            return 1
        pending.append((source_notebook, target_notebook, regenerated_text, content_changed))

    for source_notebook, target_notebook, regenerated_text, content_changed in pending:
        if not target_notebook.exists() or content_changed:
            target_notebook.parent.mkdir(parents=True, exist_ok=True)
            target_notebook.write_text(regenerated_text, encoding="utf-8")
        record_workflow_state(
            state,
            source_notebook=source_notebook,
            tracked_notebook=target_notebook,
            root=ROOT,
            operation="sync",
        )
    save_workflow_state(ROOT, state)

    print("Notebook sync complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
