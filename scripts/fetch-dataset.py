#!/usr/bin/env python3
"""Fetch a published dataset into data/, from its source, and check it.

The repository holds no data of anyone else's. What it holds, under
datasets/, is a recipe per dataset: the paper it comes with, where its files
are published, the fingerprint (sha256) each file must have, and how to turn
them into the tables the comparison scripts read. This script follows the
recipe: it downloads each file from its source into data/<name>/sources/,
refuses one whose fingerprint differs, and writes the tables to data/<name>/,
which git ignores. Whoever runs it fetches the data for their own use, under
the licence its authors gave it, which the recipe states.

    python3 scripts/fetch-dataset.py datasets/piedrafita-2020.toml
    python3 scripts/fetch-dataset.py --check datasets/*.toml

`--check` reads every recipe and checks that it says all a recipe must,
without downloading anything; the CI runs it.

Needs Python 3.11 or newer and nothing outside its standard library.
"""

import argparse
import csv
import hashlib
import io
import re
import shutil
import sys
import tomllib
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REQUIRED = ["name", "title", "authors", "source", "doi", "licence", "files", "outputs"]


def read_recipe(path):
    recipe = tomllib.loads(Path(path).read_text())
    missing = [key for key in REQUIRED if key not in recipe]
    if missing:
        raise SystemExit(f"{path}: the recipe says nothing of {', '.join(missing)}")
    for f in recipe["files"]:
        for key in ("name", "url", "sha256"):
            if key not in f:
                raise SystemExit(f"{path}: a file without `{key}`")
        if not re.fullmatch(r"[0-9a-f]{64}", f["sha256"]):
            raise SystemExit(f"{path}: {f['name']}: a sha256 is 64 hexadecimal digits")
    names = {f["name"] for f in recipe["files"]}
    for out in recipe["outputs"]:
        if out.get("convert") not in CONVERTERS:
            raise SystemExit(f"{path}: no converter named `{out.get('convert')}`; there are {sorted(CONVERTERS)}")
        if out.get("from") not in names:
            raise SystemExit(f"{path}: {out.get('name')} is made from a file the recipe does not fetch")
    return recipe


def fetch(recipe):
    directory = ROOT / "data" / recipe["name"]
    sources = directory / "sources"
    sources.mkdir(parents=True, exist_ok=True)
    for f in recipe["files"]:
        target = sources / f["name"]
        if not target.exists() or sha256(target) != f["sha256"]:
            print(f"fetching {f['url']}")
            request = urllib.request.Request(f["url"], headers={"User-Agent": "cellweave fetch-dataset"})
            with urllib.request.urlopen(request, timeout=120) as response:
                target.write_bytes(response.read())
        got = sha256(target)
        if got != f["sha256"]:
            target.unlink()
            raise SystemExit(f"{f['name']}: its fingerprint is {got}, not the {f['sha256']} the recipe expects; nothing kept")
    for out in recipe["outputs"]:
        CONVERTERS[out["convert"]](sources / out["from"], directory / out["name"])
        print(f"wrote {(directory / out['name']).relative_to(ROOT)}")
    print(f"\n{recipe['title']}\n{recipe['authors']} {recipe['source']}. https://doi.org/{recipe['doi']}\nPublished under {recipe['licence']}; cite it if you use it.")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ── Reading a spreadsheet with the standard library ──────────────────────────

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"


def sheet_rows(path, sheet):
    """Every row of one sheet of an .xlsx file, as lists of cell values from
    column A, None where a cell is empty."""
    with zipfile.ZipFile(path) as book:
        workbook = ET.fromstring(book.read("xl/workbook.xml"))
        rels = ET.fromstring(book.read("xl/_rels/workbook.xml.rels"))
        targets = {r.get("Id"): r.get("Target") or "" for r in rels}
        sheets = workbook.find("m:sheets", NS)
        rid = next((s.get(REL) for s in (sheets if sheets is not None else []) if s.get("name") == sheet), None)
        if rid is None:
            raise SystemExit(f"{path}: no sheet named {sheet}")
        target = targets[rid].lstrip("/")
        target = target if target.startswith("xl/") else "xl/" + target
        shared = []
        if "xl/sharedStrings.xml" in book.namelist():
            for si in ET.fromstring(book.read("xl/sharedStrings.xml")).findall("m:si", NS):
                shared.append("".join(t.text or "" for t in si.iter(f"{{{NS['m']}}}t")))
        data = ET.fromstring(book.read(target)).find("m:sheetData", NS)
        grid = {}
        for row in data if data is not None else []:
            for cell in row:
                ref = cell.get("r") or ""
                letters = re.match(r"[A-Z]*", ref)
                col = 0
                for ch in letters.group() if letters else "":
                    col = col * 26 + ord(ch) - 64
                kind, value = cell.get("t"), cell.find("m:v", NS)
                raw = value.text if value is not None else None
                if kind == "inlineStr":
                    text = "".join(t.text or "" for t in cell.iter(f"{{{NS['m']}}}t"))
                elif raw is None:
                    continue
                elif kind == "s":
                    text = shared[int(raw)]
                elif kind in ("str", "b", "e"):
                    text = raw
                else:
                    number = float(raw)
                    text = int(number) if number.is_integer() else number
                grid[(int(row.get("r") or 0), col)] = text
    rows = max(r for r, _ in grid)
    cols = max(c for _, c in grid)
    return [[grid.get((r, c)) for c in range(1, cols + 1)] for r in range(1, rows + 1)]


# ── Converters: from what a paper published to the table a script reads ─────


def piedrafita_2020_lrig1(source, target):
    """The sheet `Lrig1_InducedCloneSizes` of Piedrafita et al. (2020),
    Supplementary Data 10: four blocks side by side, one per time, a clone per
    row. One row per clone: days since marking (six months read as 182), mouse,
    clone, and its basal, suprabasal and terminally differentiated cells."""
    days = {"10 days": 10, "30 days": 30, "84 days": 84, "6 months": 182}
    rows = sheet_rows(source, "Lrig1_InducedCloneSizes")
    labels = [c for c in rows[4] if c is not None]
    starts = [i for i, h in enumerate(rows[6]) if h == "Mouse ID"]
    with open(target, "w", newline="") as f:
        out = csv.writer(f, lineterminator="\n")
        out.writerow(["days", "mouse", "clone", "basal", "suprabasal", "differentiated"])
        for label, s in zip(labels, starts):
            for r in rows[7:]:
                r = r + [None] * (s + 6 - len(r))
                if r[s] is None or r[s + 3] is None:
                    continue
                out.writerow([days[label], r[s], r[s + 1], int(r[s + 3]), int(r[s + 4]), int(r[s + 5])])


def copy(source, target):
    """The file as published."""
    shutil.copyfile(source, target)


CONVERTERS = {
    "piedrafita-2020-lrig1": piedrafita_2020_lrig1,
    "copy": copy,
}


def main():
    parser = argparse.ArgumentParser(description="Fetch a published dataset into data/, from its source, and check it.")
    parser.add_argument("--check", action="store_true", help="read the recipes, download nothing")
    parser.add_argument("recipes", nargs="+")
    args = parser.parse_args()
    recipes = [read_recipe(path) for path in args.recipes]
    if args.check:
        print(f"datasets: {len(recipes)} recipes, each complete")
        return
    for recipe in recipes:
        fetch(recipe)


if __name__ == "__main__":
    main()
