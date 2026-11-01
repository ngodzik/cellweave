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

With `--volume`, the runs grew in a volume: their radii are those of balls of
the same volume, and `--micrometres` says how many micrometres a voxel stands
for, so that their radii are also given in micrometres. A second table then
sets the end of the runs against the measured spheroids' limiting structure:
their radius and structure at the latest day each seeding density was imaged
at its largest, day 18 for 5,000 and 10,000 cells, day 21 for 2,500; and the
runs' over their last `--last` snapshots. A third reads, band by band of
outer radius, the radius at which Greenspan's structural relation says each
core first formed, R sqrt(1 - 3 eta^2 + 2 eta^3): with oxygen alone,
consumed at one rate, it is the same at every size.

With `--radius`, the authors' time-lapse radii as well: how fast the median
spheroid of each seeding density grows between 200 and 320 micrometres, in
micrometres a day, against how fast each run grows over the same radii, in
micrometres a step, and the hours a step that makes the two paces one.
"""

import argparse
import csv
import glob
import json
import math
import os
import statistics
from collections import defaultdict

QUIESCENCE_THRESHOLD = 0.45
BEYOND = 1.34
# The radii over which a spheroid grows at a steady pace, past its first days
# and before it slows to its end, in micrometres.
PACE_FROM, PACE_TO = 200.0, 320.0
# Bands of outer radius, in micrometres, over which the onset a core implies
# is read.
BAND = 25


def measured(path, line):
    """(R, phi, eta) of every spheroid of one cell line."""
    with open(path, newline="") as f:
        return [(float(r["R"]), float(r["ϕ"]), float(r["η"])) for r in csv.DictReader(f) if r["CellLine"] == line]


def limiting(path, line):
    """(R, phi, eta) of every spheroid imaged when its seeding was largest."""
    when = {"2500": "21", "5000": "18", "10000": "18"}
    with open(path, newline="") as f:
        return [
            (float(r["R"]), float(r["ϕ"]), float(r["η"]))
            for r in csv.DictReader(f)
            if r["CellLine"] == line and when.get(r["InitialCondition"]) == r["Day"]
        ]


def slope(points):
    """The least-squares slope of (x, y) points."""
    mx = statistics.mean(x for x, _ in points)
    my = statistics.mean(y for _, y in points)
    return sum((x - mx) * (y - my) for x, y in points) / sum((x - mx) ** 2 for x, _ in points)


def pace(path, line):
    """Per seeding density, how fast the median radius grows between PACE_FROM
    and PACE_TO, in micrometres a day."""
    radii = defaultdict(list)
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            value = float(r["R"])
            if r["CellLine"] == line and not math.isnan(value):
                radii[(int(r["InitialCondition"]), float(r["Day"]))].append(value)
    out = {}
    for seeding in sorted({s for s, _ in radii}):
        curve = sorted((day, statistics.median(v)) for (s, day), v in radii.items() if s == seeding)
        within = [(day, R) for day, R in curve if PACE_FROM <= R <= PACE_TO]
        if len(within) > 2:
            out[seeding] = slope(within)
    return out


def radius(size, volume):
    """The radius of a disc of `size` pixels, or with `volume` of a ball."""
    if volume:
        return (3 * size / (4 * math.pi)) ** (1 / 3)
    return math.sqrt(size / math.pi)


def simulated(directory, hypoxia, volume=False):
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
        R = radius(total, volume)
        inner = radius(max(total - cycling, 0), volume)
        core = radius(necrotic, volume)
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


def implied_onset(points, scale=1.0, band=None):
    """By bands of outer radius, in micrometres once multiplied by `scale`:
    how many points have a core, their median eta, and the median radius at
    which Greenspan's structural relation says their core first formed,
    Rc = R sqrt(1 - 3 eta^2 + 2 eta^3). Oxygen alone, consumed at one rate,
    keeps that radius the same at every size."""
    band = band or BAND
    out = defaultdict(list)
    for p in points:
        R, eta = p[0] * scale, p[2]
        if eta > 0 and not math.isnan(R):
            out[int(R // band) * band].append((eta, R * math.sqrt(max(0.0, 1 - 3 * eta**2 + 2 * eta**3))))
    return {b: (len(v), statistics.median(e for e, _ in v), statistics.median(r for _, r in v)) for b, v in out.items()}


def main():
    parser = argparse.ArgumentParser(description="Set spheroid runs' structure against measured spheroids'.")
    parser.add_argument("--data", required=True)
    parser.add_argument("--line", default="983b")
    parser.add_argument("--hypoxia-index", type=int, default=2, help="position of the hypoxia node among the network's nodes")
    parser.add_argument("--volume", action="store_true", help="the runs grew in a volume: radii of balls")
    parser.add_argument("--micrometres", type=float, default=1.0, help="micrometres a pixel or voxel stands for")
    parser.add_argument("--last", type=int, default=3, help="snapshots at the end of each run read as its limiting structure")
    parser.add_argument("--radius", help="the authors' time-lapse radii, to set the runs' pace against theirs")
    parser.add_argument("runs", nargs="+")
    args = parser.parse_args()

    real = structure(measured(args.data, args.line))
    per_run = [simulated(d, args.hypoxia_index, args.volume) for d in args.runs]
    model = structure([p for run in per_run for p in run])
    unit = "voxels" if args.volume else "px"
    print("| | necrosis from Rc | separated | arrest from Rq | Rq / Rc | cycling rim past 1.34 Rc | rim / Rc |")
    print("|---|---|---|---|---|---|---|")
    for name, (rc, sep, rq, rim, n) in ((f"spheroids {args.line} (um)", real), (f"runs ({unit})", model)):
        print(f"| {name} | {rc:.0f} | {sep:.0%} | {rq:.0f} | {rq / rc:.2f} | {rim:.0f} (n={n}) | {rim / rc:.2f} |")
    if args.micrometres != 1.0:
        rc, _, rq, rim, _ = model
        um = args.micrometres
        print(f"| runs (um, {um} a {unit[:-1]}) | {rc * um:.0f} | | {rq * um:.0f} | | {rim * um:.0f} | |")
    if args.radius:
        real_pace = pace(args.radius, args.line)
        # Run by run: runs that brake differently do not share a pace, and
        # pooled they would give one that none of them has.
        run_paces = []
        for run in per_run:
            steps = [(p[3], p[0] * args.micrometres) for p in run if PACE_FROM <= p[0] * args.micrometres <= PACE_TO]
            run_paces.append(slope(steps) if len({s for s, _ in steps}) > 1 else float("nan"))
        names = [os.path.basename(os.path.normpath(d)) for d in args.runs]
        print()
        print(f"| pace from {PACE_FROM:.0f} to {PACE_TO:.0f} um | um a day | " + " | ".join(f"hours a step, {n} at {p:.2f} um a step" for n, p in zip(names, run_paces)) + " |")
        print("|---|---|" + "---|" * len(run_paces))
        for seeding, per_day in real_pace.items():
            print(f"| spheroids {args.line} seeded from {seeding} | {per_day:.1f} | " + " | ".join(f"{p / per_day * 24:.2f}" for p in run_paces) + " |")
    if not args.volume:
        return
    real_onset = implied_onset(measured(args.data, args.line))
    run_onset = implied_onset([p for run in per_run for p in run], args.micrometres)
    print()
    print("| outer radius (um) | spheroids with a core | their eta | the onset it implies (um) | runs' snapshots with a core | their eta | the onset it implies (um) |")
    print("|---|---|---|---|---|---|---|")
    for b in sorted(set(real_onset) | set(run_onset)):
        cells = []
        for table in (real_onset, run_onset):
            n, eta, rc = table.get(b, (0, float("nan"), float("nan")))
            cells += [f"{n}", f"{eta:.2f}" if n else "", f"{rc:.0f}" if n else ""]
        print(f"| {b} to {b + BAND} | " + " | ".join(cells) + " |")
    ends = [p for run in per_run for p in run[-args.last:]]
    late = limiting(args.data, args.line)
    print()
    print("| at the end | R | phi | eta | cycling rim |")
    print("|---|---|---|---|---|")
    for name, pts, scale in ((f"spheroids {args.line} at their largest (um)", late, 1.0), (f"runs, last {args.last} snapshots (um)", ends, args.micrometres)):
        R = statistics.mean(p[0] for p in pts) * scale
        phi = statistics.mean(p[1] for p in pts)
        eta = statistics.mean(p[2] for p in pts)
        print(f"| {name} | {R:.0f} | {phi:.2f} | {eta:.2f} | {R * (1 - phi):.0f} |")


if __name__ == "__main__":
    main()
