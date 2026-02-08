import argparse
from pathlib import Path


def _write_text_exact(path: Path, content: str) -> None:
    if content and not content.endswith("\n"):
        content += "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Copy mirrored OSF wiki pages from docs/osf_wiki into the main docs/*.md files used by mkdocs.yml nav.",
    )
    parser.add_argument(
        "--docs-dir",
        default="docs",
        help="Docs directory (default: docs).",
    )
    parser.add_argument(
        "--osf-wiki-dir",
        default="docs/osf_wiki",
        help="Mirrored OSF wiki directory (default: docs/osf_wiki).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would change without writing files.",
    )
    args = parser.parse_args()

    docs_dir = Path(args.docs_dir)
    osf_dir = Path(args.osf_wiki_dir)
    if not docs_dir.exists():
        raise SystemExit(f"docs dir not found: {docs_dir}")
    if not osf_dir.exists():
        raise SystemExit(
            f"OSF wiki dir not found: {osf_dir}. Run tools/mirror_osf_wiki.py first."
        )

    mapping = {
        # Home
        "home.md": "index.md",
        # Installation
        "1. Installation - Python.md": "Installation-Python.md",
        "2. Installation - R.md": "Installation-R.md",
        # Usage
        "3. Plotting Substitutions.md": "Plotting-Substitutions.md",
        "4. Plotting Indels.md": "Plotting-Indels.md",
        "5. Plotting Dinucleotides.md": "Plotting-Dinucleotides.md",
        "6. Plotting a Sample Portrait.md": "Plotting-a-Sample-Portrait.md",
        "7. Example Program.md": "Example-Program.md",
    }

    missing_sources = []
    for src_name in mapping:
        src_path = osf_dir / src_name
        if not src_path.exists():
            missing_sources.append(src_name)
    if missing_sources:
        raise SystemExit(
            "Missing expected OSF wiki pages:\n"
            + "\n".join(f"- {name}" for name in missing_sources)
        )

    changed = 0
    for src_name, dst_name in mapping.items():
        src_path = osf_dir / src_name
        dst_path = docs_dir / dst_name
        content = src_path.read_text(encoding="utf-8")

        if args.dry_run:
            print(f"{src_path} -> {dst_path}")
            continue

        _write_text_exact(dst_path, content)
        changed += 1

    if not args.dry_run:
        print(f"Updated {changed} docs pages from {osf_dir}")
        print("Preview with: python3 -m mkdocs serve")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
