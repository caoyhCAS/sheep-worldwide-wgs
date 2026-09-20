"""Reject duplicate/missing/overlapping sample assignments instead of silently dropping them."""
import argparse
import csv
import gzip
import json
from pathlib import Path


def read_ids(path, columns=1):
    ids = []
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != columns or (columns == 2 and fields[0] != fields[1]):
                raise ValueError(f"{path}:{number}: expected {'identical FID IID' if columns == 2 else 'one sample ID'}")
            ids.append(fields[-1])
    if not ids or len(ids) != len(set(ids)):
        raise ValueError(f"{path}: empty list or duplicate sample IDs")
    return set(ids)


def vcf_samples(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM\t"):
                ids = line.rstrip("\n").split("\t")[9:]
                if not ids or len(ids) != len(set(ids)):
                    raise ValueError("VCF has missing or duplicate sample IDs")
                return set(ids)
    raise ValueError("VCF has no #CHROM header")


def validate(vcf, populations, unrelated=None, check_ld=False):
    available = vcf_samples(vcf)
    allowed = read_ids(unrelated, 2) if unrelated else available
    if not allowed <= available:
        raise ValueError(f"Unrelated samples absent from VCF: {sorted(allowed - available)}")
    assigned = set()
    report = {"vcf_samples": len(available), "allowed_samples": len(allowed), "populations": {}}
    with open(populations, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"population", "keep"} <= set(reader.fieldnames or []):
            raise ValueError("Population manifest requires population and keep columns")
        for row in reader:
            name = row["population"]
            if name in report["populations"]:
                raise ValueError(f"Duplicate population: {name}")
            members = read_ids(row["keep"])
            if not members <= allowed:
                raise ValueError(f"Population {name}: samples absent from allowed cohort: {sorted(members - allowed)}")
            if members & assigned:
                raise ValueError(f"Population assignments overlap: {sorted(members & assigned)}")
            assigned |= members
            report["populations"][name] = {"n": len(members)}
            if check_ld:
                selected = read_ids(row["ld_keep"])
                if len(selected) != 3 or not selected <= members:
                    raise ValueError(f"Population {name}: LD list must contain exactly three member samples")
                report["populations"][name]["ld_samples"] = sorted(selected)
    if not report["populations"]:
        raise ValueError("Empty population manifest")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--populations", required=True)
    parser.add_argument("--unrelated")
    parser.add_argument("--check-ld", action="store_true")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = validate(args.vcf, args.populations, args.unrelated, args.check_ld)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
