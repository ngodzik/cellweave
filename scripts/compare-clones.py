#!/usr/bin/env python3
"""Set the clones of traced runs against a lineage tracing experiment's.

For every time the experiment counted clones, take the marked clones of each
run at that time after marking, keep those with a cell in the layer, as the
experiment keeps those with a basal cell, and compare: how many cells the
surviving clones hold on average, the shares holding more than one, two and
three times that mean, and the largest gap between the two distributions of
sizes (the Kolmogorov-Smirnov distance, 0 for the same distribution, 1 for
two that do not overlap).

    python3 scripts/compare-clones.py --data data/piedrafita-2020/lrig1-clones.csv \\
        --at 1000 --hours-per-step 0.6 output-1 output-2 ...

`--at` is the step the runs marked their cells at. `--hours-per-step` says
what a step stands for; `--division-days` instead makes each run's own
division interval, read from its `measurements.csv` after marking, stand for
that many days, which is how chapter 12 calibrates a sheet. Each run
directory holds the `traced.csv` and `measurements.csv` the run wrote. The rows of the experiment are pooled over its mice, and those of the
runs over the runs.
"""

import argparse
import csv
import os
import statistics
import sys


def experiment(path):
    """Basal cells of every clone, by days since marking."""
    sizes = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            sizes.setdefault(int(row["days"]), []).append(int(row["basal"]))
    return sizes


def steps_per_division(directory, at):
    """How often a living cell divided after marking, in steps: the living,
    averaged over the steps from marking on, times their number, over the
    births in them, as a run's summary reckons it."""
    rows = []
    with open(os.path.join(directory, "measurements.csv"), newline="") as f:
        for row in csv.DictReader(f):
            if int(row["mcs"]) >= at:
                rows.append(row)
    living = statistics.mean(int(r["cycling"]) + int(r["quiescent"]) for r in rows)
    births = int(rows[-1]["born"]) - int(rows[0]["born"])
    return living * (int(rows[-1]["mcs"]) - int(rows[0]["mcs"])) / births


def progenitor_share(directory, at):
    """The share of cells in the layer that have not committed to
    differentiate, averaged over the steps from marking on."""
    shares = []
    with open(os.path.join(directory, "measurements.csv"), newline="") as f:
        for row in csv.DictReader(f):
            if int(row["mcs"]) >= at:
                living = int(row["cycling"]) + int(row["quiescent"])
                layer = living + int(row["differentiating"])
                if layer:
                    shares.append(living / layer)
    return statistics.mean(shares)


def runs(dirs, at, hours, division_days, days):
    """Cells of every marked clone still in the layer, by days since marking,
    pooled over the runs, each read at the step written down nearest to the
    day. A run's step stands for `hours`, or, given `division_days`, for
    whatever makes its own division interval that many days."""
    sizes = {d: [] for d in days}
    for directory in dirs:
        written = {}
        with open(os.path.join(directory, "traced.csv"), newline="") as f:
            for row in csv.DictReader(f):
                written.setdefault(int(row["mcs"]), []).append(int(row["cells"]))
        h = hours if division_days is None else division_days * 24 / steps_per_division(directory, at)
        for d in days:
            want = at + d * 24 / h
            nearest = min(written, key=lambda mcs: abs(mcs - want))
            if abs(nearest - want) > 24 / h:
                print(f"{directory}: nothing written within a day of {d} days (step {want:.0f})", file=sys.stderr)
                continue
            sizes[d].extend(written[nearest])
    return sizes


def ks(a, b):
    """The Kolmogorov-Smirnov distance between two samples of integers."""
    a, b = sorted(a), sorted(b)
    gap, i, j = 0.0, 0, 0
    for x in sorted(set(a) | set(b)):
        while i < len(a) and a[i] <= x:
            i += 1
        while j < len(b) and b[j] <= x:
            j += 1
        gap = max(gap, abs(i / len(a) - j / len(b)))
    return gap


def shares(sizes):
    """Shares of clones larger than one, two and three times the mean."""
    mean = statistics.mean(sizes)
    return [sum(s > k * mean for s in sizes) / len(sizes) for k in (1, 2, 3)]


def main():
    parser = argparse.ArgumentParser(description="Set traced runs' clones against a lineage tracing experiment's.")
    parser.add_argument("--data", required=True)
    parser.add_argument("--at", type=int, required=True)
    when = parser.add_mutually_exclusive_group(required=True)
    when.add_argument("--hours-per-step", type=float)
    when.add_argument("--division-days", type=float)
    parser.add_argument(
        "--progenitors",
        type=float,
        help="the experiment's share of basal cells that are progenitors, to set beside the runs'",
    )
    parser.add_argument("runs", nargs="+")
    args = parser.parse_args()

    measured = experiment(args.data)
    simulated = runs(args.runs, args.at, args.hours_per_step, args.division_days, sorted(measured))
    print("| days | clones, mouse | in the runs | mean cells, mouse | in the runs | ratio | above 1, 2, 3 times the mean, mouse | in the runs | KS |")
    print("|---|---|---|---|---|---|---|---|---|")
    for days in sorted(measured):
        mouse = [s for s in measured[days] if s >= 1]
        model = simulated[days]
        if not model:
            print(f"| {days} | {len(mouse)} | 0 | {statistics.mean(mouse):.2f} | | | | | |")
            continue
        m, s = statistics.mean(mouse), statistics.mean(model)
        fm = ", ".join(f"{x:.2f}" for x in shares(mouse))
        fs = ", ".join(f"{x:.2f}" for x in shares(model))
        print(f"| {days} | {len(mouse)} | {len(model)} | {m:.2f} | {s:.2f} | {s / m:.2f} | {fm} | {fs} | {ks(mouse, model):.2f} |")

    if args.progenitors is not None:
        model = statistics.mean(progenitor_share(d, args.at) for d in args.runs)
        print(f"\nShare of the layer not committed to differentiate: {args.progenitors:.2f} in the experiment, {model:.2f} in the runs.")


if __name__ == "__main__":
    main()
