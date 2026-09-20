#!/usr/bin/env python3
"""Auditable adapters for selected msab353 methods, not the authors' original code.

The default action only prints a JSON plan. External programs are executed only
with --execute, in a newly created run directory, without a command shell.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone


MODULES = ("psmc", "smcpp", "fastsimcoal", "qpdstat", "abba", "popgen", "smartpca")
PROGRAMS = {
    "psmc": "psmc", "smcpp": "smc++", "fastsimcoal": "fsc26",
    "qpdstat": "qpDstat", "abba": "python3", "popgen": "python3",
    "smartpca": "smartpca",
}
REPORTED_VERSIONS = {
    "psmc": "not reported", "smcpp": "1.14.0.dev0", "fastsimcoal": "2.6",
    "qpdstat": "ADMIXTOOLS 7.0.1", "abba": "not reported",
    "popgen": "not reported", "smartpca": "EIGENSOFT 6.0.1",
}
SOURCES = {
    "psmc": "https://github.com/lh3/psmc/blob/master/README",
    "smcpp": "https://github.com/popgenmethods/smcpp",
    "fastsimcoal": "https://cmpg.unibe.ch/software/fastsimcoal26/man/fastsimcoal26.pdf",
    "qpdstat": "https://github.com/DReichLab/AdmixTools/blob/master/README.Dstatistics",
    "abba": "https://github.com/simonhmartin/genomics_general/blob/master/ABBABABAwindows.py",
    "popgen": "https://github.com/simonhmartin/genomics_general/blob/master/popgenWindows.py",
    "smartpca": "https://github.com/DReichLab/EIG/blob/master/POPGEN/README",
}
CONFIG_FIELDS = {
    "psmc": {"psmcfa"},
    "smcpp": {"smc_files"},
    "fastsimcoal": {"template", "estimation", "observed_sfs", "seed"},
    "qpdstat": {"parameter_file", "chromosome_26_support_confirmed", "supported_build_evidence"},
    "smartpca": {"parameter_file"},
    "abba": {"script", "genotypes", "populations_file", "genotype_format", "populations", "step_bp", "threads", "ploidy_file"},
    "popgen": {"script", "genotypes", "populations_file", "genotype_format", "populations", "step_bp", "threads", "ploidy_file"},
}


def text_value(value, label):
    if not isinstance(value, str) or not value or any(c in value for c in "\x00\r\n"):
        raise ValueError(f"{label} must be a nonempty single-line string")
    return value


def required(config, key):
    if key not in config:
        raise ValueError(f"required configuration field is missing: {key}")
    return config[key]


def positive_integer(value, label, upper=None):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{label} must be a positive integer")
    if upper is not None and value > upper:
        raise ValueError(f"{label} must be at most {upper}")
    return value


def absolute_path(value, base, label):
    path = Path(text_value(value, label)).expanduser()
    return str((path if path.is_absolute() else base / path).resolve())


def parse_parameter_file(path, module):
    """Restrict path-bearing settings so an external tool cannot overwrite inputs.

    These intentionally small allowlists are not complete third-party parsers.
    Unrecognized keys are rejected rather than silently changing their meaning.
    """
    input_keys = {"genotypename", "snpname", "indivname"}
    if module == "qpdstat":
        input_keys.add("popfilename")
        scalar_keys = {"f4mode", "printsd", "blgsize", "numchrom", "inbreed"}
        output_keys = {}
    else:
        scalar_keys = {
            "numoutevec", "numoutlieriter", "numoutlierevec", "outliersigmathresh",
            "numchrom", "usenorm", "altnormstyle", "missingmode", "lsqproject",
        }
        output_keys = {"evecoutname": "results.evec", "evaloutname": "results.eval"}
    values = {}
    dependencies = {}
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        key, separator, value = line.partition(":")
        key, value = key.strip().lower(), value.strip()
        if not separator or not value:
            raise ValueError(f"{path}:{line_number}: expected key: value")
        if key in values:
            raise ValueError(f"{path}:{line_number}: duplicate parameter {key}")
        if key not in input_keys | scalar_keys | output_keys.keys():
            raise ValueError(f"{path}:{line_number}: unsupported parameter {key}; review manually")
        if key in input_keys:
            dependencies[key] = absolute_path(value, path.parent, key)
            # Keep spaces and shell metacharacters out of the legacy parameter
            # parser by linking the read-only input to a fixed local filename.
            value = f"input.{key}"
        elif key in output_keys:
            value = output_keys[key]
        elif not re.fullmatch(r"(?:YES|NO|[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)", value, re.I):
            raise ValueError(f"{key}: expected a numeric or YES/NO scalar")
        values[key] = value
    missing = sorted(input_keys - values.keys())
    if missing:
        raise ValueError(f"{path}: missing parameter-file inputs: {', '.join(missing)}")
    if module == "qpdstat" and values.get("f4mode", "NO").upper() != "NO":
        raise ValueError("qpdstat adapter requires f4mode: NO (D statistics, not f4)")
    if values.get("numchrom") != "26":
        raise ValueError(f"{module} requires explicit numchrom: 26 for sheep autosomal inputs; the human default is unsuitable")
    values.update(output_keys)
    content = "".join(f"{key}: {value}\n" for key, value in values.items())
    return content, dependencies


def build_plan(module, config, run_dir, base_dir=None):
    """Read configuration/available parameter files; never create files or run tools."""
    if module not in MODULES:
        raise ValueError(f"unknown module: {module}")
    if not isinstance(config, dict):
        raise ValueError("module configuration must be a JSON object")
    unknown = set(config) - CONFIG_FIELDS[module] - {"executable", "requested_version"}
    if unknown:
        raise ValueError(f"unsupported {module} configuration fields: {', '.join(sorted(unknown))}")
    base = Path(base_dir or Path.cwd()).resolve()
    directory = absolute_path(str(run_dir), base, "run_dir")
    program = text_value(config.get("executable", PROGRAMS[module]), "executable")
    if os.sep in program:
        program = absolute_path(program, base, "executable")
    plan = {
        "schema_version": 1,
        "status": "plan_only",
        "method_status": "method-based reconstruction; not original author code",
        "module": module,
        "reported_version": REPORTED_VERSIONS[module],
        "requested_version": config.get("requested_version", REPORTED_VERSIONS[module]),
        "version_verification": "not performed; install and verify the requested tool separately",
        "cli_reference": SOURCES[module],
        "cwd": directory,
        "argv": [program],
        "required_files": [],
        "staged_files": [],
        "linked_files": [],
        "generated_files": {},
        "expected_outputs": [],
        "stdout": "stdout.log",
        "stderr": "stderr.log",
        "blockers": [],
        "notes": ["All run outputs stay in a new directory; existing directories are refused."],
    }

    def input_path(key):
        path = absolute_path(required(config, key), base, key)
        add_input(path)
        return path

    def add_input(path):
        if not any(item["path"] == path for item in plan["required_files"]):
            plan["required_files"].append({"path": path, "exists": Path(path).is_file()})

    def output(path, kind="file"):
        plan["expected_outputs"].append({"path": path, "kind": kind})

    if module == "psmc":
        data = input_path("psmcfa")
        plan["argv"] += ["-N25", "-t15", "-r5", "-p", "4+25*2+4+6", "-o", "result.psmc", data]
        output("result.psmc")
        plan["notes"] += [
            "Input must be a callable-mask-aware diploid psmcfa, not a SNP-only VCF.",
            "Raw output is scaled PSMC units, NOT years. Plotting requires mutation rate 1.51e-8 and generation time 3 years.",
            "The reported highest-depth-three-individuals selection is an upstream responsibility.",
        ]
    elif module == "smcpp":
        data_files = required(config, "smc_files")
        if not isinstance(data_files, list) or not data_files:
            raise ValueError("smc_files must be a nonempty list of prepared SMC++ input files")
        paths = [absolute_path(p, base, "smc_files entry") for p in data_files]
        for path in paths:
            add_input(path)
        plan["argv"] += ["estimate", "-o", "estimate", "1.51e-8", *paths]
        output("estimate/model.final.json")
        plan["notes"] += [
            "Input preparation, distinguished individuals, accessible/mappability masks and sample choice are not inferred.",
            "The model is NOT a years-scaled plot. SMC++ plotting/scaling is a separate step with generation time 3 years.",
        ]
    elif module == "fastsimcoal":
        seed = positive_integer(required(config, "seed"), "seed", 1000000)
        template, estimation, observed = input_path("template"), input_path("estimation"), input_path("observed_sfs")
        plan["staged_files"] = [
            {"source": template, "target": "model.tpl"},
            {"source": estimation, "target": "model.est"},
            {"source": observed, "target": "model_DSFS.obs"},
        ]
        plan["argv"] += [
            "-t", "model.tpl", "-e", "model.est", "-n", "1000000", "-d", "-M",
            "-l", "25", "-L", "65", "--multiSFS", "--seed", str(seed),
        ]
        output("model", "directory")
        if Path(estimation).is_file():
            active = "\n".join(line.split("//", 1)[0] for line in Path(estimation).read_text(encoding="utf-8").splitlines())
            if not re.search(r"\breference\b", active, re.I):
                plan["blockers"].append("-l 25 requires a reference parameter in the supplied .est (fsc26 manual); none found")
        plan["notes"] += [
            "Runs ONE optimization, not the paper's 100 restarts. Repeat with recorded distinct seeds and new run directories.",
            "Supply a validated unfolded multidimensional SFS including monomorphic sites; staged filename is model_DSFS.obs.",
            "Specify mutation rate 1.51e-8 in the user-reviewed demographic model; this wrapper does not infer or edit model biology.",
            "The three demographic scenarios, search ranges, 2647-block bootstrap and 100 x 20 bootstrap optimizations are not generated here.",
        ]
    elif module in {"qpdstat", "smartpca"}:
        par = input_path("parameter_file")
        if Path(par).is_file():
            content, dependencies = parse_parameter_file(Path(par), module)
            plan["generated_files"]["analysis.par"] = content
            for key, path in dependencies.items():
                add_input(path)
                plan["linked_files"].append({"source": path, "target": f"input.{key}"})
        else:
            plan["blockers"].append("parameter file is unavailable; its input dependencies cannot yet be checked")
        plan["argv"] += ["-p", "analysis.par"]
        if module == "qpdstat":
            plan["stdout"] = "results.txt"
            output("results.txt")
            confirmed = config.get("chromosome_26_support_confirmed", False)
            if not isinstance(confirmed, bool):
                raise ValueError("chromosome_26_support_confirmed must be a JSON boolean")
            if not confirmed:
                plan["blockers"].append("Stock ADMIXTOOLS v7.0.1 ignores numchrom and drops sheep chr23-26; explicitly verify a patched/portable build before setting chromosome_26_support_confirmed=true")
            else:
                evidence = text_value(required(config, "supported_build_evidence"), "supported_build_evidence")
                plan["supported_build_evidence"] = evidence
                requested_version = text_value(required(config, "requested_version"), "requested_version")
                if re.fullmatch(r"(?:ADMIXTOOLS\s+)?v?7\.0\.1", requested_version, re.I):
                    plan["blockers"].append("requested_version still identifies unmodified v7.0.1; explicitly identify the audited patched or portable build")
            plan["notes"].append("Explicit W X Y Z tests must be in popfilename; |Z| > 3 is the reported significance criterion.")
            plan["notes"].append("The supplied parameter file must explicitly set numchrom: 26 and EIGENSTRAT inputs must use sheep autosomes 1-26.")
            plan["notes"].append("Chromosome support is caller-attested, not automatically verified. A patched/portable build is a documented departure from the paper's reported stock v7.0.1.")
        else:
            output("results.evec")
            output("results.eval")
            plan["notes"].append("Input EIGENSTRAT conversion and sample pruning are upstream. The supplied parameter file must explicitly set numchrom: 26 and inputs must use sheep autosomes 1-26.")
        plan["notes"].append("A restricted, rewritten analysis.par uses local links to original inputs; unsupported keys fail closed. The original parameter file is unchanged.")
    else:
        script, genotypes, populations = input_path("script"), input_path("genotypes"), input_path("populations_file")
        add_input(str(Path(script).with_name("genomics.py")))
        step = required(config, "step_bp")
        if isinstance(step, bool) or not isinstance(step, int) or step not in (20000, 50000):
            raise ValueError("step_bp must explicitly select 20000 (reported command) or 50000 (reported prose)")
        genotype_format = required(config, "genotype_format")
        if genotype_format not in ("phased", "pairs", "haplo", "diplo"):
            raise ValueError("genotype_format must be phased, pairs, haplo or diplo; VCF is not accepted directly")
        threads = positive_integer(config.get("threads", 1), "threads")
        plan["argv"] += [
            script, "-g", genotypes, "-o", "windows.csv", "-f", genotype_format,
            "--popsFile", populations, "-w", "100000", "-m", "500", "-s", str(step), "-T", str(threads),
        ]
        names = required(config, "populations")
        expected_count = 4 if module == "abba" else 2
        if not isinstance(names, list) or len(names) != expected_count:
            raise ValueError(f"populations must contain exactly {expected_count} names")
        names = [text_value(name, "population name") for name in names]
        if len(set(names)) != len(names) or any(name.startswith("-") or any(c.isspace() for c in name) for name in names):
            raise ValueError("population names must be distinct, without whitespace or a leading '-'")
        flags = ["-P1", "-P2", "-P3", "-O"] if module == "abba" else ["-p", "-p"]
        for flag, name in zip(flags, names):
            plan["argv"] += [flag, name]
        if module == "popgen":
            plan["argv"] += ["--analysis", "popDist", "popPairDist"]
        if "ploidy_file" in config:
            plan["argv"] += ["--ploidyFile", input_path("ploidy_file")]
        output("windows.csv")
        plan["notes"] += [
            "Conflict retained: paper prose gives a 50 kb step; the printed command gives -s 20000. The chosen step is explicit.",
            "Supply the full genomics_general checkout (including genomics.py) and record its commit; no third-party code is vendored.",
            "Genotype preparation, ploidy and callable-site denominator require review. Missing/inaccessible sites must not be treated as invariant reference.",
        ]
        if module == "abba":
            plan["notes"].append("populations order is P1 reference (Menz), P2 recipient, P3 donor, O outgroup. Ancestral-state pseudo-outgroups must be supplied, never invented.")
        else:
            plan["notes"].append("Run twice for donor/recipient and donor/reference. For absolute pi/dxy, provide callable invariant sites as well as variants.")
    for item in plan["required_files"]:
        if not item["exists"]:
            plan["blockers"].append(f"missing input file: {item['path']}")
    if Path(directory).exists():
        plan["blockers"].append(f"run directory already exists: {directory}")
    return plan


def execute_plan(plan):
    """Execute a locally constructed plan once; keep failed runs for inspection."""
    if plan["blockers"]:
        raise ValueError("execution blocked: " + "; ".join(plan["blockers"]))
    for item in plan["required_files"]:
        if not Path(item["path"]).is_file():
            raise ValueError(f"input disappeared: {item['path']}")
    executable = shutil.which(plan["argv"][0])
    if executable is None:
        raise ValueError(f"executable not found or not executable: {plan['argv'][0]}")
    directory = Path(plan["cwd"])
    # mkdir(exist_ok=False) also refuses an existing empty directory/symlink.
    directory.mkdir(parents=True, exist_ok=False)
    record = dict(plan)
    record["status"] = "running"
    record["resolved_executable"] = executable
    record["started_utc"] = datetime.now(timezone.utc).isoformat()
    record_path = directory / "run.json"

    def save_record():
        record_path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    save_record()
    try:
        for item in plan["staged_files"]:
            shutil.copyfile(item["source"], directory / item["target"])
        for item in plan["linked_files"]:
            (directory / item["target"]).symlink_to(item["source"])
        for name, content in plan["generated_files"].items():
            (directory / name).write_text(content, encoding="utf-8")
        argv = [executable, *plan["argv"][1:]]
        with (directory / plan["stdout"]).open("x", encoding="utf-8") as out, (directory / plan["stderr"]).open("x", encoding="utf-8") as err:
            result = subprocess.run(argv, cwd=directory, stdin=subprocess.DEVNULL, stdout=out, stderr=err, shell=False, check=False)
        record["returncode"] = result.returncode
        absent = []
        for item in plan["expected_outputs"]:
            path = directory / item["path"]
            if not (path.is_dir() if item["kind"] == "directory" else path.is_file()):
                absent.append(item["path"])
        record["missing_expected_outputs"] = absent
        record["status"] = "completed" if result.returncode == 0 and not absent else "failed"
        record["output_validation"] = "file/directory existence only; scientific results not validated"
    except Exception as exc:
        record["status"] = "failed"
        record["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        record["finished_utc"] = datetime.now(timezone.utc).isoformat()
        record["files_present"] = sorted(str(p.relative_to(directory)) for p in directory.rglob("*") if p.is_file())
        save_record()
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="JSON config; file paths are relative to the invocation working directory")
    parser.add_argument("--module", choices=MODULES, required=True)
    parser.add_argument("--run-dir", required=True, help="fresh output directory; must not already exist")
    parser.add_argument("--execute", action="store_true", help="run the external tool (default prints a non-mutating JSON plan)")
    args = parser.parse_args(argv)
    try:
        with open(args.config, encoding="utf-8") as handle:
            configuration = json.load(handle)
        plan = build_plan(args.module, required(configuration, args.module), args.run_dir)
        result = execute_plan(plan) if args.execute else plan
        print(json.dumps(result, indent=2))
        return 0 if result["status"] in ("plan_only", "completed") else 1
    except (OSError, ValueError, TypeError) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    sys.exit(main())
