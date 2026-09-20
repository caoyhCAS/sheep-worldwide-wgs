"""Pass one -V per input to GATK without shell interpolation."""
import argparse
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("vcfs", nargs="+")
    args = parser.parse_args()
    command = ["gatk", "--java-options", "-Xmx24g", "CombineGVCFs", "-R", args.reference, "-O", args.output]
    for path in args.vcfs:
        command.extend(["-V", path])
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
