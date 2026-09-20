"""Run FastQC in an isolated temporary directory, preserving arbitrary FASTQ names."""
import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--stem", required=True)
    parser.add_argument("--html", required=True)
    parser.add_argument("--zip", required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="msab353-fastqc-") as directory:
        source = Path(args.input).resolve()
        linked = Path(directory) / (args.stem + (".fastq.gz" if source.name.endswith(".gz") else ".fastq"))
        linked.symlink_to(source)
        subprocess.run(["fastqc", "--outdir", directory, "--noextract", str(linked)], check=True)
        for extension, destination in [("html", args.html), ("zip", args.zip)]:
            Path(destination).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(Path(directory) / f"{args.stem}_fastqc.{extension}", destination)


if __name__ == "__main__":
    main()
