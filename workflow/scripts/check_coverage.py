"""Whole-reference depth summary and strict > threshold admission gate.

Zero-depth positions must be included by `samtools depth -aa`. This denominator
and samtools defaults are reconstruction choices, not specified by the paper.
"""
import argparse
import json
import math
import sys
from pathlib import Path


def measure(lines, expected_bases):
    count = total = covered = 0
    for line in lines:
        fields = line.rstrip().split("\t")
        if len(fields) != 3:
            raise ValueError("Expected exactly three samtools depth columns for one BAM")
        depth = int(fields[2])
        if depth < 0:
            raise ValueError("Negative depth")
        count += 1
        total += depth
        covered += depth > 0
    if count != expected_bases or count == 0:
        raise ValueError(f"Expected {expected_bases} reference positions including zeros; saw {count}")
    return {"positions": count, "covered_positions": covered, "coverage_fraction": covered / count,
            "mean_depth": total / count, "definition": "samtools depth -aa, all FASTA bases, tool defaults"}


def gate(summary, threshold):
    depth = float(summary["mean_depth"])
    if not math.isfinite(depth) or not depth > threshold:
        raise ValueError(f"Mean depth {depth} does not satisfy strict > {threshold}; revise cohort explicitly")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="command", required=True)
    measure_parser = subs.add_parser("measure")
    measure_parser.add_argument("--sample", required=True)
    measure_parser.add_argument("--fai", required=True)
    measure_parser.add_argument("--output", required=True)
    gate_parser = subs.add_parser("gate")
    gate_parser.add_argument("--input", required=True)
    gate_parser.add_argument("--threshold", type=float, required=True)
    gate_parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.command == "measure":
        with open(args.fai, encoding="utf-8") as handle:
            expected = sum(int(line.split("\t")[1]) for line in handle if line.strip())
        report = measure(sys.stdin, expected)
        report["sample"] = args.sample
        result = json.dumps(report, indent=2) + "\n"
    else:
        with open(args.input, encoding="utf-8") as handle:
            gate(json.load(handle), args.threshold)
        result = f"mean_depth > {args.threshold}\n"
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(result, encoding="utf-8")


if __name__ == "__main__":
    main()
