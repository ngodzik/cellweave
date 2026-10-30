#!/usr/bin/env python3
"""How fast clones grow when every cell that leaves is replaced by a neighbour.

A voter model on a hexagonal lattice, wrapped round at its edges: at each
update a cell picked at random is lost, and the daughter of one of its six
neighbours, picked at random, takes its place. Every cell starts as a clone of
its own. Prints the clones still there and their mean size after a number of
divisions per cell, the time unit chapter 24 compares in.

It stands for the most local replacement a sheet could have, with no
coordination between sisters, as a yardstick for the sheet's clones.

    python3 scripts/voter-model.py
"""

import collections
import random

SIDE = 100
STEPS = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]


def main():
    rng = random.Random(1)
    cells = SIDE * SIDE
    clone = list(range(cells))
    for divisions in range(1, 76):
        for _ in range(cells):
            i = rng.randrange(cells)
            dx, dy = STEPS[rng.randrange(6)]
            x, y = i % SIDE, i // SIDE
            clone[i] = clone[(x + dx) % SIDE + ((y + dy) % SIDE) * SIDE]
        if divisions in (4, 12, 35, 75):
            alive = len(collections.Counter(clone))
            print(f"after {divisions} divisions per cell: {alive} clones, {cells / alive:.1f} cells each")


if __name__ == "__main__":
    main()
