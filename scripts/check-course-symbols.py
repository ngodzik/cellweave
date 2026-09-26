#!/usr/bin/env python3
"""Every code identifier the course names must still exist in the code.

Scans docs/*.md for backticked tokens that look like identifiers, a name with
an underscore, a double colon, or a leading capital, and greps each in the Rust
sources. A chapter that names a symbol which was renamed or removed fails the
check, which is the most common way a course drifts from its program.
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
SOURCES = [ROOT / "crates", ROOT / "src"]

TOKEN = re.compile(r"`([^`\n]+)`")
IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_:]*$")


def looks_like_code(token: str) -> bool:
    if not IDENT.match(token):
        return False
    return "_" in token or "::" in token or token[0].isupper()


def exists(symbol: str) -> bool:
    needle = symbol.split("::")[-1]
    result = subprocess.run(
        ["grep", "-rqw", "--include=*.rs", needle, *map(str, SOURCES)],
        check=False,
    )
    return result.returncode == 0


def main() -> int:
    missing: list[tuple[Path, str]] = []
    seen: set[str] = set()
    for page in sorted(DOCS.glob("*.md")):
        for token in TOKEN.findall(page.read_text()):
            if not looks_like_code(token) or token in seen:
                continue
            seen.add(token)
            if not exists(token):
                missing.append((page.relative_to(ROOT), token))
    if missing:
        print("the course names symbols the code no longer has:")
        for page, token in missing:
            print(f"  {page}: `{token}`")
        return 1
    print(f"course symbols: {len(seen)} checked, all present")
    return 0


if __name__ == "__main__":
    sys.exit(main())
