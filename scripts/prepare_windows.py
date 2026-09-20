#!/usr/bin/env python3
"""Join native VCFtools windows without assuming that row orders coincide.

VCFtools BIN_START/BIN_END are retained as 1-based inclusive coordinates.
New reconstruction glue, not an original study script.
"""
import argparse
import csv
import math
import sys
from pathlib import Path


def read_windows(path, value_column):
    result = {}
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"CHROM", "BIN_START", "BIN_END", value_column}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: requires columns {sorted(required)}")
        for line, row in enumerate(reader, 2):
            chrom = row["CHROM"].strip()
            start, end = int(row["BIN_START"]), int(row["BIN_END"])
            if not chrom or start < 1 or end < start:
                raise ValueError(f"{path}:{line}: invalid 1-based inclusive interval")
            key = (chrom, start, end)
            if key in result:
                raise ValueError(f"{path}:{line}: duplicate interval {key}")
            value = row[value_column].strip()
            if value.lower() in {"", ".", "na", "nan", "-nan"}:
                value = "NA"
            elif not math.isfinite(float(value)):
                raise ValueError(f"{path}:{line}: infinite value")
            result[key] = value
    if not result:
        raise ValueError(f"{path}: contains no windows")
    return result


def join_windows(tables, mode="exact"):
    keys = [set(table) for table in tables]
    common = set.intersection(*keys)
    if mode == "exact" and any(keyset != keys[0] for keyset in keys[1:]):
        raise ValueError("Window sets differ; investigate or explicitly use --join intersection")
    if not common:
        raise ValueError("No common windows")
    # Preserve the first file's chromosome ordering, not lexical chr1/chr10 order.
    return [key for key in tables[0] if key in common]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    rod = sub.add_parser("rod", help="Join VCFtools .windowed.pi tables")
    rod.add_argument("--landrace", required=True)
    rod.add_argument("--improved", required=True)
    pbs = sub.add_parser("pbs", help="Join three VCFtools .windowed.weir.fst tables")
    pbs.add_argument("--cha-efr", required=True)
    pbs.add_argument("--cha-hus", required=True)
    pbs.add_argument("--efr-hus", required=True)
    pbs.add_argument("--fst-column", choices=["MEAN_FST", "WEIGHTED_FST"], required=True,
                     help="Explicit estimator choice; not the paper's original R FST estimator")
    for child in [rod, pbs]:
        child.add_argument("--join", choices=["exact", "intersection"], default="exact")
        child.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    if args.command == "rod":
        paths, column = [args.landrace, args.improved], "PI"
        names = ["pi_landrace", "pi_improved"]
    else:
        paths, column = [args.cha_efr, args.cha_hus, args.efr_hus], args.fst_column
        names = ["fst_cha_efr", "fst_cha_hus", "fst_efr_hus"]
    try:
        if Path(args.output).resolve() in {Path(p).resolve() for p in paths}:
            raise ValueError("Output must not overwrite an input")
        tables = [read_windows(path, column) for path in paths]
        keys = join_windows(tables, args.join)
        dest = Path(args.output)
        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["chrom", "start", "end", *names])
            for key in keys:
                writer.writerow([*key, *[table[key] for table in tables]])
        print(f"Wrote {len(keys)} windows; input counts={[len(t) for t in tables]}; join={args.join}", file=sys.stderr)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
