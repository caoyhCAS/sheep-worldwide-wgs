#!/usr/bin/env python3
"""Auditable, standard-library statistics reconstructed from msab353 methods.

These are not the authors' original scripts. See window_stats.md for assumptions.
Input tables use normalized, unique (chrom, start, end) interval keys.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import statistics
import sys
from typing import Iterable


KEY_COLUMNS = ("chrom", "start", "end")
MISSING = {"", ".", "na", "nan", "+nan", "-nan"}


def number(value: str, field: str, *, missing: bool = False) -> float | None:
    """Parse a finite number; recognized missing tokens require explicit consent."""
    if value.strip().lower() in MISSING:
        if missing:
            return None
        raise ValueError(f"{field}: missing value is not allowed")
    try:
        result = float(value)
    except ValueError as error:
        raise ValueError(f"{field}: not a numeric value: {value!r}") from error
    if not math.isfinite(result):
        raise ValueError(f"{field}: infinity/nonfinite value is not allowed")
    return result


def read_windows(path: str | Path, required: Iterable[str]) -> tuple[list[str], list[dict[str, str]]]:
    """Validate normalized TSV keys, including duplicates and malformed row widths."""
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError("duplicate column names")
        absent = set(KEY_COLUMNS).union(required).difference(fields)
        if absent:
            raise ValueError(f"missing columns: {', '.join(sorted(absent))}")
        rows = []
        seen = set()
        for line, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"line {line}: wrong number of TSV fields")
            chrom = row["chrom"]
            if not chrom or chrom != chrom.strip():
                raise ValueError(f"line {line}: chromosome must be nonempty with no surrounding whitespace")
            try:
                start, end = int(row["start"]), int(row["end"])
            except ValueError as error:
                raise ValueError(f"line {line}: interval coordinates must be integers") from error
            if start < 1 or end < start:
                raise ValueError(f"line {line}: expected 1 <= start <= end (1-based inclusive)")
            key = (chrom, start, end)
            if key in seen:
                raise ValueError(f"line {line}: duplicate window {key}")
            seen.add(key)
            rows.append(row)
    if not rows:
        raise ValueError("input table contains no windows")
    return fields, rows


def write_table(path: str | Path, fields: list[str], rows: Iterable[dict]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "NA" if value is None else value for key, value in row.items()})


def append_fields(fields: list[str], new: list[str]) -> list[str]:
    """Avoid silently replacing source columns when reprocessing an output file."""
    collision = set(fields).intersection(new)
    if collision:
        raise ValueError(f"output columns already exist: {', '.join(sorted(collision))}")
    return fields + new


def pbs(fst_cha_efr: float, fst_cha_hus: float, fst_efr_hus: float, *, negative_fst: str = "error") -> float:
    if negative_fst not in {"error", "zero"}:
        raise ValueError("negative_fst must be error or zero")
    transformed = []
    for value in (fst_cha_efr, fst_cha_hus, fst_efr_hus):
        if not math.isfinite(value) or value >= 1:
            raise ValueError("PBS requires finite FST < 1; FST=1 has infinite drift length")
        if value < 0:
            if negative_fst == "error":
                raise ValueError("negative FST requires explicit --negative-fst zero")
            value = 0.0
        transformed.append(-math.log1p(-value))
    return (transformed[0] + transformed[1] - transformed[2]) / 2


def rod(pi_landrace: float, pi_improved: float) -> float | None:
    if not all(math.isfinite(value) and value >= 0 for value in (pi_landrace, pi_improved)):
        raise ValueError("nucleotide diversity must be finite and nonnegative")
    if pi_landrace == 0:
        return None
    result = 1 - pi_improved / pi_landrace
    if not math.isfinite(result):
        raise ValueError("ROD overflow: check diversity units and magnitudes")
    return result


def benjamini_hochberg(values: list[float | None]) -> list[float | None]:
    finite = []
    for index, value in enumerate(values):
        if value is None:
            continue
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError("BH input must be probabilities in [0,1] or None")
        finite.append((value, index))
    finite.sort()
    adjusted: list[float | None] = [None] * len(values)
    running = 1.0
    n = len(finite)
    for reverse_index in range(n - 1, -1, -1):
        value, original_index = finite[reverse_index]
        running = min(running, value * n / (reverse_index + 1))
        adjusted[original_index] = running
    return adjusted


def clean_fd(fd: float | None, d: float | None) -> tuple[float | None, str]:
    if any(value is not None and not math.isfinite(value) for value in (fd, d)):
        raise ValueError("fd and D must be finite or missing")
    if d is not None and not -1 <= d <= 1:
        raise ValueError("D must be in [-1,1] or missing")
    if fd is None or d is None:
        return None, "missing"
    if fd < 0 or d < 0 or fd > 1:
        return 0.0, "set_to_zero"
    return fd, "retained"


def fd_z_bh(values: list[float | None]) -> tuple[list[float | None], list[float | None], list[float | None], str]:
    """Upper-tail normal approximation, sample SD, BH family = finite windows."""
    finite = [value for value in values if value is not None]
    if any(not math.isfinite(value) for value in finite):
        raise ValueError("fd values must be finite or None")
    if len(finite) < 2:
        absent = [None] * len(values)
        return absent.copy(), absent.copy(), absent.copy(), "insufficient_windows"
    mean = statistics.mean(finite)
    sd = statistics.stdev(finite)
    if sd == 0:
        absent = [None] * len(values)
        return absent.copy(), absent.copy(), absent.copy(), "zero_variance"
    zscores = [None if value is None else (value - mean) / sd for value in values]
    pvalues = [None if z is None else 0.5 * math.erfc(z / math.sqrt(2)) for z in zscores]
    return zscores, pvalues, benjamini_hochberg(pvalues), "ok"


def empirical_threshold(values: list[float], top_fraction: float) -> float:
    """R type 7 / NumPy linear quantile; selection is >= threshold, including ties."""
    if not math.isfinite(top_fraction) or not 0 < top_fraction < 1:
        raise ValueError("top fraction must lie strictly between 0 and 1")
    if not values or any(not math.isfinite(value) for value in values):
        raise ValueError("quantile requires at least one finite value")
    ordered = sorted(values)
    position = (len(ordered) - 1) * (1 - top_fraction)
    lower = math.floor(position)
    upper = math.ceil(position)
    fraction = position - lower
    # Weighted interpolation avoids overflow from (upper_value - lower_value).
    return (1 - fraction) * ordered[lower] + fraction * ordered[upper]


def ils_probability(tract_bp: float, divergence_years: float, generation_years: float,
                    recombination_rate: float, branch_factor: float = 2.0) -> dict[str, float | str]:
    """Gamma(shape=2, rate=r*t) survival, keeping log probability under underflow."""
    parameters = (tract_bp, divergence_years, generation_years, recombination_rate, branch_factor)
    if not all(math.isfinite(value) for value in parameters):
        raise ValueError("ILS parameters must be finite")
    if tract_bp < 0 or min(divergence_years, generation_years, recombination_rate, branch_factor) <= 0:
        raise ValueError("tract length must be nonnegative; other ILS parameters must be positive")
    branch_generations = branch_factor * (divergence_years / generation_years)
    rate = recombination_rate * branch_generations
    x = tract_bp * rate
    if not math.isfinite(rate) or rate <= 0 or not math.isfinite(x):
        raise ValueError("ILS parameter product overflow/underflow; rescale inputs")
    expected_length = 1 / rate
    if not math.isfinite(expected_length):
        raise ValueError("ILS expected tract length overflow; check input magnitudes")
    log_probability = math.log1p(x) - x
    probability = math.exp(log_probability)
    return {
        "tract_bp": tract_bp,
        "divergence_years": divergence_years,
        "generation_years": generation_years,
        "recombination_rate_per_bp_per_generation": recombination_rate,
        "branch_factor": branch_factor,
        "branch_generations": branch_generations,
        "expected_length_bp": expected_length,
        "log_probability": log_probability,
        "probability": probability,
        "probability_status": "underflow" if probability == 0 else "finite",
    }


def run(args: argparse.Namespace) -> None:
    if args.command == "ils":
        result = ils_probability(args.tract_bp, args.divergence_years, args.generation_years,
                                 args.recombination_rate, args.branch_factor)
        write_table(args.output, list(result), [result])
        return

    required = {
        "pbs": ["fst_cha_efr", "fst_cha_hus", "fst_efr_hus"],
        "rod": ["pi_landrace", "pi_improved"],
        "fd": ["fd", "D"],
        "outliers": [getattr(args, "column", "score")],
    }[args.command]
    fields, rows = read_windows(args.input, required)
    if args.command == "pbs":
        fields = append_fields(fields, ["pbs", "negative_fst_policy"])
        for row in rows:
            row["pbs"] = pbs(*(number(row[key], key) for key in required), negative_fst=args.negative_fst)
            row["negative_fst_policy"] = args.negative_fst
    elif args.command == "rod":
        fields = append_fields(fields, ["rod", "rod_status"])
        for row in rows:
            row["rod"] = rod(*(number(row[key], key) for key in required))
            row["rod_status"] = "zero_denominator" if row["rod"] is None else "ok"
    elif args.command == "fd":
        fields = append_fields(fields, ["fd_clean", "clean_status", "z", "p_upper", "q_bh", "test_status", "test_n", "tail", "sd"])
        for row in rows:
            row["fd_clean"], row["clean_status"] = clean_fd(number(row["fd"], "fd", missing=True),
                                                          number(row["D"], "D", missing=True))
        values = [row["fd_clean"] for row in rows]
        zscores, pvalues, qvalues, status = fd_z_bh(values)
        test_n = sum(value is not None for value in values)
        for row, z, p, q in zip(rows, zscores, pvalues, qvalues):
            row.update(z=z, p_upper=p, q_bh=q, test_status="missing" if row["fd_clean"] is None else status,
                       test_n=test_n, tail=args.tail, sd=args.sd)
    elif args.command == "outliers":
        fields = append_fields(fields, ["quantile_threshold", "is_outlier", "top_fraction", "quantile_method", "finite_n"])
        values = [number(row[args.column], args.column, missing=args.missing == "exclude") for row in rows]
        finite = [value for value in values if value is not None]
        threshold = empirical_threshold(finite, args.top_fraction)
        for row, value in zip(rows, values):
            row.update(quantile_threshold=threshold, is_outlier=None if value is None else int(value >= threshold),
                       top_fraction=args.top_fraction, quantile_method="linear_type7_ge", finite_n=len(finite))
    write_table(args.output, fields, rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    pbs_parser = commands.add_parser("pbs", help="PBS from one exactly aligned table of three pairwise FSTs")
    pbs_parser.add_argument("--negative-fst", choices=["error", "zero"], default="error")
    commands.add_parser("rod", help="1 - pi_improved/pi_landrace; zero denominator produces NA")
    fd_parser = commands.add_parser("fd", help="Clean fd, reconstructed window Z transform and BH")
    fd_parser.add_argument("--tail", choices=["upper"], required=True,
                           help="Explicitly accept reconstructed upper-tail normal approximation")
    fd_parser.add_argument("--sd", choices=["ddof1"], required=True,
                           help="Explicitly accept reconstructed sample-standard-deviation choice")
    outliers = commands.add_parser("outliers", help="Empirical upper quantile, including all boundary ties")
    outliers.add_argument("--column", required=True)
    outliers.add_argument("--top-fraction", required=True, type=float)
    outliers.add_argument("--missing", choices=["error", "exclude"], default="error")
    for command in (pbs_parser, commands.choices["rod"], fd_parser, outliers):
        command.add_argument("--input", required=True)
        command.add_argument("--output", required=True)
    ils = commands.add_parser("ils", help="ILS gamma-shape-2 survival with log probability")
    ils.add_argument("--tract-bp", required=True, type=float)
    ils.add_argument("--divergence-years", required=True, type=float)
    ils.add_argument("--generation-years", required=True, type=float)
    ils.add_argument("--recombination-rate", required=True, type=float)
    ils.add_argument("--branch-factor", type=float, default=2.0,
                     help="Total branch length multiplier, default 2 (two lineages); verify biological model")
    ils.add_argument("--output", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if getattr(args, "input", None) and Path(args.input).resolve() == Path(args.output).resolve():
            raise ValueError("input and output paths must differ")
        run(args)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())
