#!/usr/bin/env python3
"""Reproducible candidate sample selection, not historical study sample lists.

Required TSV: sample, population, mean_depth. IDs must be unique.
For PSMC select highest coverage; for LD select a seeded random subset.
Relatedness must be resolved upstream; this helper does not compute kinship.
"""
import argparse
import csv
import math
import random
import sys
from pathlib import Path


def select(rows, mode, number, seed=None):
    if number < 1:
        raise ValueError("number must be positive")
    if mode not in {"highest-depth", "random"}:
        raise ValueError("Unknown selection mode")
    if mode == "random" and seed is None:
        raise ValueError("Random selection requires a recorded seed")
    groups, seen = {}, set()
    for row in rows:
        sample, pop = row["sample"].strip(), row["population"].strip()
        depth = float(row["mean_depth"])
        if not sample or not pop or sample in seen:
            raise ValueError("Empty sample/population or duplicate sample ID")
        if not math.isfinite(depth) or depth < 0:
            raise ValueError("mean_depth must be finite and nonnegative")
        seen.add(sample)
        groups.setdefault(pop, []).append({**row, "sample": sample, "population": pop})
    if not groups:
        raise ValueError("Empty sample table")
    chosen = []
    for pop in sorted(groups):
        group = sorted(groups[pop], key=lambda row: row["sample"])
        if len(group) < number:
            raise ValueError(f"{pop}: requires {number} individuals, only {len(group)} available")
        if mode == "highest-depth":
            chosen.extend(sorted(group, key=lambda row: (-float(row["mean_depth"]), row["sample"]))[:number])
        else:
            # Per-population generator prevents unrelated new groups changing a selection.
            rng = random.Random(f"{seed}:{pop}")
            chosen.extend(sorted(rng.sample(group, number), key=lambda row: row["sample"]))
    return chosen


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--mode", choices=["highest-depth", "random"], required=True)
    parser.add_argument("--number", type=int, default=3)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    try:
        if Path(args.input).resolve() == Path(args.output).resolve():
            raise ValueError("Output must not overwrite input")
        with open(args.input, encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"sample", "population", "mean_depth"}.issubset(reader.fieldnames or []):
                raise ValueError("Required columns: sample, population, mean_depth")
            chosen = select(list(reader), args.mode, args.number, args.seed)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w", encoding="utf-8", newline="") as handle:
            fields = ["sample", "population", "mean_depth", "selection_mode", "seed"]
            writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            for row in chosen:
                writer.writerow({**row, "selection_mode": args.mode, "seed": args.seed if args.seed is not None else "NA"})
    except (ValueError, OSError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
