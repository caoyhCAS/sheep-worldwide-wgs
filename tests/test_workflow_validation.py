"""Validate reconstruction gates and format every real rule without running tools.

Empty placeholders in the DAG test are never interpreted as biological inputs:
Snakemake dry-run checks only path dependencies and command formatting.
"""
import importlib.util
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[1]


def helper(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "workflow" / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_coverage_counts_zero_positions_and_strict_gate():
    coverage = helper("check_coverage")
    summary = coverage.measure(["1\t1\t0\n", "1\t2\t20\n"], 2)
    assert summary["mean_depth"] == 10
    assert summary["coverage_fraction"] == 0.5
    with pytest.raises(ValueError, match="strict"):
        coverage.gate(summary, 10)
    coverage.gate({"mean_depth": 10.01}, 10)
    with pytest.raises(ValueError):
        coverage.gate({"mean_depth": float("nan")}, 10)
    with pytest.raises(ValueError, match="Expected 3"):
        coverage.measure(["1\t1\t5\n"], 3)


def test_cohort_validator_rejects_unknown_overlap_and_wrong_ld(tmp_path):
    cohort = helper("validate_cohort")
    vcf = tmp_path / "input.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\n")
    a = tmp_path / "a.txt"
    b = tmp_path / "b.txt"
    a.write_text("A\nB\nC\n")
    b.write_text("D\n")
    manifest = tmp_path / "pops.tsv"
    manifest.write_text(f"population\tkeep\tld_keep\nP1\t{a}\t{a}\nP2\t{b}\t{b}\n")
    assert cohort.validate(vcf, manifest)["populations"]["P1"]["n"] == 3
    with pytest.raises(ValueError, match="exactly three"):
        cohort.validate(vcf, manifest, check_ld=True)
    b.write_text("A\n")
    with pytest.raises(ValueError, match="overlap"):
        cohort.validate(vcf, manifest)
    b.write_text("ABSENT\n")
    with pytest.raises(ValueError, match="absent"):
        cohort.validate(vcf, manifest)


def real_config(tmp_path):
    def placeholder(name, content=""):
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return str(path)

    placeholder("config/config.test.yaml", "modules: [demo]\ndemo:\n  targets: [results/demo/provenance.json]\n")
    reference = placeholder("reference.fa", ">1\nA\n")
    for suffix in [".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"]:
        placeholder("reference.fa" + suffix, "1\t1\t3\t1\t2\n" if suffix == ".fai" else "")
    dictionary = placeholder("reference.dict")
    r1 = placeholder("sample_R1.fq.gz")
    r2 = placeholder("sample_R2.fq.gz")
    samples = placeholder("samples.tsv", f"sample\tr1\tr2\tmean_depth\nS1\t{r1}\t{r2}\t20\n")
    vcf = placeholder("input.vcf.gz")
    placeholder("input.vcf.gz.tbi")
    rename = placeholder("contigs.tsv", "1\t1\n")
    a = placeholder("a.txt", "A\nB\nC\n")
    b = placeholder("b.txt", "D\nE\nF\n")
    pops = placeholder("populations.tsv", f"population\tkeep\tld_keep\nP1\t{a}\t{a}\nP2\t{b}\t{b}\n")
    unrelated = placeholder("unrelated.keep", "A A\nB B\nC C\nD D\nE E\nF F\n")
    return {
        "modules": ["wgs", "population", "selection"], "threads": 2,
        "reference": {"fasta": reference, "dict": dictionary},
        "wgs": {"samples": samples, "trimmomatic_steps": ["MINLEN:1"], "missing_annotation_policy": "pass"},
        "analysis": {"vcf": vcf, "contig_map": rename, "autosomes": ["1"], "populations": pops},
        "population": {"unrelated_keep": unrelated, "roh": {"enabled": True, "resolution_note": "Syntax-only test choice, not scientific recommendation", "resolved_args": ["--homozyg-kb", "200"]}, "ld": {"enabled": True}},
        "selection": {"fst_window_bp": 50000, "fst_step_bp": 10000, "fst_step_resolution": "Use Methods for syntax-only test", "global_fst": True, "comparisons": {"P1_vs_P2": ["P1", "P2"]}},
    }


def snakemake_command():
    executable = Path(sys.executable).parent / "snakemake"
    if executable.is_file():
        return str(executable)
    found = shutil.which("snakemake")
    if not found:
        pytest.skip("Snakemake is required for DAG validation; install requirements-dev.txt")
    return found


def test_every_real_rule_dry_run_formats_shell(tmp_path):
    config = real_config(tmp_path)
    path = tmp_path / "full.yaml"
    path.write_text(yaml.safe_dump(config))
    result = subprocess.run([snakemake_command(), "--snakefile", str(ROOT / "workflow" / "Snakefile"), "--configfile", str(path), "--cores", "2", "--dry-run", "--printshellcmds"], cwd=tmp_path, text=True, capture_output=True)
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    expected = ["raw_fastqc", "trim_reads", "align_reads", "remove_duplicates", "measure_coverage", "coverage_gate", "haplotype_caller", "combine_gvcfs", "genotype_gvcfs", "select_variant_type", "hard_filter", "pass_variants", "prepare_autosomal_vcf", "validate_analysis_cohort", "plink_import", "plink_unrelated", "population_qc", "ld_prune_list", "apply_ld_pruning", "ibs_distance", "population_heterozygosity", "population_pi", "population_roh", "population_ld", "selection_fst", "selection_global_fst", "selection_pi"]
    for rule in expected:
        assert f"rule {rule}:" in output, rule
    assert "rule individual_pi:" in output
    assert "--fst-window-size 50000 --fst-window-step 10000" in output
    assert "--REMOVE_DUPLICATES true" in output
    assert "--restrict-alleles-to BIALLELIC" in output
    assert "--distance square 1-ibs" in output
    assert "-MaxDist 300" in output
    assert "--filter-name QUAL30 --filter-expression 'QUAL < 30.0'" in output
    assert "--filter-name MQRankSum --filter-expression 'MQRankSum < -12.5'" in output
    assert "--filter-name FS200 --filter-expression 'FS > 200.0'" in output
    assert "||" not in output


def test_unresolved_selection_step_is_blocked(tmp_path):
    config = real_config(tmp_path)
    config["selection"]["fst_step_bp"] = None
    path = tmp_path / "unresolved.yaml"
    path.write_text(yaml.safe_dump(config))
    result = subprocess.run([snakemake_command(), "--snakefile", str(ROOT / "workflow" / "Snakefile"), "--configfile", str(path), "--dry-run"], cwd=tmp_path, text=True, capture_output=True)
    assert result.returncode != 0
    assert "Required configuration: fst_step_bp" in result.stdout + result.stderr
