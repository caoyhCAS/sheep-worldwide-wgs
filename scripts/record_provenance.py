#!/usr/bin/env python3
"""Hash explicitly named inputs/configs; do not scan environment or credentials."""
import argparse
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True, help="Repeat for config/input files")
    parser.add_argument("--label", default="method-based reconstruction; not original code")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if Path(args.output).resolve() in {Path(p).resolve() for p in args.input}:
        parser.error("Output cannot overwrite an input")
    manifest = {"study_doi": "10.1093/molbev/msab353", "label": args.label,
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "python": platform.python_version(),
                "inputs": [{"path": p, "bytes": Path(p).stat().st_size, "sha256": sha256(p)} for p in args.input]}
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
