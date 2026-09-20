#!/usr/bin/env python3
"""Assert synthetic demonstration outputs. Never asserts biological reproduction."""
import csv
import json
import math
from pathlib import Path


def table(name):
    with Path(f"results/demo/{name}.tsv").open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main():
    ils = table("ils")[0]
    assert math.isclose(float(ils["probability"]), 0.0008650084490313926, rel_tol=1e-12)
    assert ils["probability_status"] == "finite"
    rod = table("rod")
    assert float(rod[0]["rod"]) == 0.75
    assert rod[2]["rod"] == "NA"
    fd = table("fd")
    assert float(fd[1]["fd_clean"]) == 0
    assert fd[-1]["q_bh"] == "NA"
    assert [row["is_outlier"] for row in table("pbs_outliers")] == ["0", "0", "0", "1"]
    provenance = json.loads(Path("results/demo/provenance.json").read_text())
    assert len(provenance["inputs"]) == 6
    print("Synthetic demo outputs verified; this is not full study reproduction.")


if __name__ == "__main__":
    main()
