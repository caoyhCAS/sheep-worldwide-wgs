"""Small standard-library validators shared by the Snakefile and tests."""
import csv
import re


def required(mapping, dotted):
    current = mapping
    for key in dotted.split("."):
        if not isinstance(current, dict) or key not in current or current[key] in (None, ""):
            raise ValueError(f"Required configuration: {dotted}")
        current = current[key]
    return current


def safe_name(value, label):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise ValueError(f"{label} must contain only letters, digits, '_', '.' or '-' and start alphanumeric: {value!r}")
    return value


def positive_integer(value, label):
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return value


def manifest(path, columns, key):
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not set(columns).issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: required columns {columns}")
        rows = {}
        for row in reader:
            if any(not row.get(c) for c in columns):
                raise ValueError(f"{path}: empty required cell")
            if row[key] in rows:
                raise ValueError(f"{path}: duplicate {key} {row[key]}")
            rows[row[key]] = row
    if not rows:
        raise ValueError(f"Empty manifest: {path}")
    return rows
