#!/usr/bin/env python3
"""Set a spheroid run's structure against measured spheroids'.

A spheroid's structure is three radii: the outer radius R, the radius inside
which cells are arrested or dead, and the radius of the necrotic core. Browning
and colleagues (eLife, 2021) measured all three on hundreds of melanoma
spheroids, and found that the structure follows the outer radius rather than
time. This script reads the snapshots of runs in a bath, works out the same
three radii from them, as radii of discs of the same area, and sets both
against each other in three numbers each:

- the radius at which a necrotic core appears, Rc: the cut between two
  observed radii that best separates spheroids with a core from those without;
- the radius at which arrest appears, Rq, fitted by least squares as Browning
  and colleagues write it before a core forms, phi = sqrt(max(0, 1 - Rq^2/R^2)),
  phi being the arrested or dead radius over R, zero below Rq;
- the rim of cycling cells, R less the arrested or dead radius, once R is past
  1.34 Rc, the measured spheroids' 320 micrometres.

    python3 scripts/compare-spheroids.py --data data/browning-2021/confocal.csv \\
        --line 983b output-spheroid-1 output-spheroid-2 ...

A run's cell is arrested when its hypoxia node is at or above the threshold
the simulation uses, 0.45, the node being named by `--hypoxia-index` in the
order of the network's nodes; cells of kind Necrotic are dead. A run's radii
are in pixels and a measured spheroid's in micrometres, a pixel standing for
about a micrometre (chapter 1). The two are compared by their ratios to their
own Rc, so that a disc and a sphere, whose cores appear at radii a factor
sqrt(2/3) apart for the same oxygen, are set side by side fairly.
"""

import argparse
import csv
import glob
import json
import math
import os
import statistics

QUIESCENCE_THRESHOLD = 0.45
BEYOND = 1.34


def measured(path, line):
    """(R, phi, eta) of every spheroid of one cell line."""
    with open(path, newline="") as f:
        return [(float(r["R"]), float(r["ϕ"]), float(r["η"])) for r in csv.DictReader(f) if r["CellLine"] == line]


def simulated(directory, hypoxia):
    """(R, phi, eta, step) from every snapshot of a run, once it has a cell."""
    out = []
    for path in sorted(glob.glob(os.path.join(directory, "snapshot_*.json"))):
        with open(path) as f:
            snap = json.load(f)
        total = cycling = necrotic = 0
        for cell in snap["cells"]:
            area = cell["volume"]
            total += area
            if cell["kind"] == "Necrotic":
                necrotic += area
            elif cell["kind"] == "Tumor" and cell["proteins"]["values"][hypoxia] < QUIESCENCE_THRESHOLD:
                cycling += area
        if total == 0:
            continue
        R = math.sqrt(total / math.pi)
        inner = math.sqrt(max(total - cycling, 0) / math.pi)
        core = math.sqrt(necrotic / math.pi)
        out.append((R, inner / R, core / R, snap["mcs"]))
    return out


def structure(points):
    """Rc, Rq and the cycling rim past 1.34 Rc, from (R, phi, eta, ...)."""
    radii = sorted({p[0] for p in points})
    cuts = [(a + b) / 2 for a, b in zip(radii, radii[1:])]
    rc = max(cuts, key=lambda c: sum((p[0] >= c) == (p[2] > 0) for p in points))
    separated = sum((p[0] >= rc) == (p[2] > 0) for p in points) / len(points)
    # phi = sqrt(max(0, 1 - Rq^2 / R^2)) before a core forms, zero below Rq:
    # the Rq that fits it best by least squares, searched to a tenth.
    before = [p for p in points if p[2] == 0]
    def misfit(rq):
        return sum((p[1] - math.sqrt(max(0.0, 1 - rq * rq / (p[0] * p[0])))) ** 2 for p in before)
    rq = min((x / 10 for x in range(0, int(rc * 10) + 1)), key=misfit)
    rims = [p[0] * (1 - p[1]) for p in points if p[0] > BEYOND * rc]
    rim = statistics.median(rims) if rims else float("nan")
    return rc, separated, rq, rim, len(rims)


def main():
    parser = argparse.ArgumentParser(description="Set spheroid runs' structure against measured spheroids'.")
    parser.add_argument("--data", required=True)
    parser.add_argument("--line", default="983b")
    parser.add_argument("--hypoxia-index", type=int, default=2, help="position of the hypoxia node among the network's nodes")
    parser.add_argument("runs", nargs="+")
    args = parser.parse_args()

    real = structure(measured(args.data, args.line))
    runs = [p for d in args.runs for p in simulated(d, args.hypoxia_index)]
    model = structure(runs)
    print("| | necrosis from Rc | separated | arrest from Rq | Rq / Rc | cycling rim past 1.34 Rc | rim / Rc |")
    print("|---|---|---|---|---|---|---|")
    for name, (rc, sep, rq, rim, n) in ((f"spheroids {args.line} (um)", real), ("runs (px)", model)):
        print(f"| {name} | {rc:.0f} | {sep:.0%} | {rq:.0f} | {rq / rc:.2f} | {rim:.0f} (n={n}) | {rim / rc:.2f} |")


if __name__ == "__main__":
    main()
