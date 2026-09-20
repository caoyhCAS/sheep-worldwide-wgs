"""Adapter validation only: these tests do not run any biological software."""

import json
from pathlib import Path
import subprocess
import sys

import pytest

from legacy_tools import build_plan, execute_plan, main


def write_file(path, content="fixture\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return str(path)


def genomic_config(tmp_path, module="abba"):
    script = write_file(tmp_path / "external tools" / "tool.py")
    write_file(tmp_path / "external tools" / "genomics.py")
    return {
        "script": script,
        "genotypes": write_file(tmp_path / "data with spaces" / "input.geno"),
        "populations_file": write_file(tmp_path / "pops.txt"),
        "genotype_format": "diplo",
        "populations": ["Menz", "recipient", "donor", "outgroup"] if module == "abba" else ["donor", "recipient"],
        "step_bp": 20000,
    }


def parameter_config(tmp_path, module="qpdstat", suffix=""):
    keys = ["genotypename", "snpname", "indivname"]
    if module == "qpdstat":
        keys.append("popfilename")
    lines = [f"{key}: {write_file(tmp_path / 'input with spaces' / key)}" for key in keys]
    lines += ["numchrom: 26"]
    if module == "smartpca":
        lines += ["evecoutname: /existing/should-not-overwrite", "evaloutname: ../outside.eval"]
    path = write_file(tmp_path / "params.par", "\n".join(lines) + "\n" + suffix)
    return {"parameter_file": path}


def test_psmc_paper_parameters_and_nonmutating_plan(tmp_path):
    data = write_file(tmp_path / "sample.psmcfa")
    plan = build_plan("psmc", {"psmcfa": data}, tmp_path / "new run")
    assert plan["argv"] == ["psmc", "-N25", "-t15", "-r5", "-p", "4+25*2+4+6", "-o", "result.psmc", data]
    assert plan["blockers"] == []
    assert not (tmp_path / "new run").exists()
    assert "not performed" in plan["version_verification"]


def test_missing_inputs_are_reported_but_never_created(tmp_path):
    plan = build_plan("psmc", {"psmcfa": "absent.psmcfa"}, tmp_path / "run", tmp_path)
    assert plan["required_files"] == [{"path": str(tmp_path / "absent.psmcfa"), "exists": False}]
    with pytest.raises(ValueError, match="execution blocked"):
        execute_plan(plan)
    assert not (tmp_path / "run").exists()


def test_unknown_options_fail_instead_of_being_silently_ignored(tmp_path):
    with pytest.raises(ValueError, match="unsupported psmc configuration"):
        build_plan("psmc", {"psmcfa": "input", "extra_args": ["-o", "/outside"]}, tmp_path / "run")


def test_smcpp_mutation_parameter_and_no_implied_year_scaling(tmp_path):
    data = write_file(tmp_path / "population.smc.gz")
    plan = build_plan("smcpp", {"smc_files": [data]}, tmp_path / "run")
    assert plan["argv"] == ["smc++", "estimate", "-o", "estimate", "1.51e-8", data]
    assert any("NOT a years-scaled" in note for note in plan["notes"])
    with pytest.raises(ValueError, match="nonempty list"):
        build_plan("smcpp", {"smc_files": []}, tmp_path / "run")


def test_fsc_one_seeded_run_and_safe_staging(tmp_path):
    cfg = {
        "template": write_file(tmp_path / "scenario I.tpl"),
        "estimation": write_file(tmp_path / "scenario I.est", "[PARAMETERS]\n1 NPOP unif 10 100000 output reference\n"),
        "observed_sfs": write_file(tmp_path / "scenario I_DSFS.obs"),
        "seed": 234,
    }
    plan = build_plan("fastsimcoal", cfg, tmp_path / "run")
    assert plan["blockers"] == []
    assert plan["argv"][-6:] == ["25", "-L", "65", "--multiSFS", "--seed", "234"]
    assert "1000000" in plan["argv"]
    assert "--runs" not in plan["argv"]
    assert [item["target"] for item in plan["staged_files"]] == ["model.tpl", "model.est", "model_DSFS.obs"]
    Path(cfg["estimation"]).write_text("// reference comment must not count\n[PARAMETERS]\n", encoding="utf-8")
    plan = build_plan("fastsimcoal", cfg, tmp_path / "run")
    assert any("reference parameter" in blocker for blocker in plan["blockers"])
    with pytest.raises(ValueError, match="seed"):
        build_plan("fastsimcoal", {**cfg, "seed": 1000001}, tmp_path / "run")


@pytest.mark.parametrize("step", [20000, 50000])
def test_introgression_step_conflict_is_user_resolved(tmp_path, step):
    cfg = genomic_config(tmp_path)
    cfg["step_bp"] = step
    plan = build_plan("abba", cfg, tmp_path / "run")
    argv = plan["argv"]
    assert argv[argv.index("-s") + 1] == str(step)
    assert argv[argv.index("-w") + 1] == "100000"
    assert argv[argv.index("-m") + 1] == "500"
    assert argv[-8:] == ["-P1", "Menz", "-P2", "recipient", "-P3", "donor", "-O", "outgroup"]


def test_no_silent_introgression_defaults_or_option_injection(tmp_path):
    cfg = genomic_config(tmp_path)
    del cfg["step_bp"]
    with pytest.raises(ValueError, match="step_bp"):
        build_plan("abba", cfg, tmp_path / "run")
    cfg["step_bp"] = 10000
    with pytest.raises(ValueError, match="20000"):
        build_plan("abba", cfg, tmp_path / "run")
    cfg["step_bp"] = 20000
    cfg["populations"][0] = "--outFile"
    with pytest.raises(ValueError, match="leading"):
        build_plan("abba", cfg, tmp_path / "run")


def test_popgen_requests_within_and_between_population_statistics(tmp_path):
    cfg = genomic_config(tmp_path, "popgen")
    cfg["ploidy_file"] = write_file(tmp_path / "ploidy.txt")
    plan = build_plan("popgen", cfg, tmp_path / "run")
    assert plan["argv"][-9:] == ["-p", "donor", "-p", "recipient", "--analysis", "popDist", "popPairDist", "--ploidyFile", cfg["ploidy_file"]]
    assert any(item["path"].endswith("genomics.py") for item in plan["required_files"])


def test_qpdstat_checks_nested_inputs_and_preserves_original_par(tmp_path):
    cfg = parameter_config(tmp_path)
    original = Path(cfg["parameter_file"]).read_text(encoding="utf-8")
    plan = build_plan("qpdstat", cfg, tmp_path / "run")
    assert len(plan["required_files"]) == 5
    assert len(plan["linked_files"]) == 4
    assert plan["generated_files"]["analysis.par"].startswith("genotypename: input.genotypename\n")
    assert plan["stdout"] == "results.txt"
    assert Path(cfg["parameter_file"]).read_text(encoding="utf-8") == original
    Path(plan["linked_files"][0]["source"]).unlink()
    plan = build_plan("qpdstat", cfg, tmp_path / "run")
    assert any("missing input file" in item for item in plan["blockers"])


def test_qpdstat_stock_version_is_blocked_even_with_numchrom26(tmp_path):
    cfg = parameter_config(tmp_path)
    plan = build_plan("qpdstat", cfg, tmp_path / "run")
    assert any("Stock ADMIXTOOLS v7.0.1" in item for item in plan["blockers"])
    with pytest.raises(ValueError, match="execution blocked"):
        execute_plan(plan)
    cfg.update({
        "chromosome_26_support_confirmed": True,
        "supported_build_evidence": "Synthetic assertion for unit test; not an audited scientific binary",
        "requested_version": "ADMIXTOOLS 7.0.1",
    })
    plan = build_plan("qpdstat", cfg, tmp_path / "run")
    assert any("unmodified v7.0.1" in item for item in plan["blockers"])
    cfg["requested_version"] = "v7.0.1 with locally reviewed numchrom patch (test fixture)"
    plan = build_plan("qpdstat", cfg, tmp_path / "run")
    assert plan["blockers"] == []
    assert "Synthetic assertion" in plan["supported_build_evidence"]


def test_qpdstat_confirmation_requires_recorded_evidence(tmp_path):
    cfg = parameter_config(tmp_path)
    cfg["chromosome_26_support_confirmed"] = True
    with pytest.raises(ValueError, match="supported_build_evidence"):
        build_plan("qpdstat", cfg, tmp_path / "run")


@pytest.mark.parametrize("suffix", ["f4mode: YES\n", "output: /tmp/uncontrolled\n", "include: other.par\n", "printsd: $(touch nope)\n"])
def test_qpdstat_rejects_unsupported_semantics(tmp_path, suffix):
    cfg = parameter_config(tmp_path, suffix=suffix)
    with pytest.raises(ValueError):
        build_plan("qpdstat", cfg, tmp_path / "run")


def test_smartpca_outputs_are_confined_to_fresh_directory(tmp_path):
    cfg = parameter_config(tmp_path, "smartpca")
    plan = build_plan("smartpca", cfg, tmp_path / "run")
    par = plan["generated_files"]["analysis.par"]
    assert "evecoutname: results.evec\n" in par
    assert "evaloutname: results.eval\n" in par
    assert "/existing/should-not-overwrite" not in par


@pytest.mark.parametrize("chromosome_line", ["", "numchrom: 22\n", "numchrom: 27\n"])
@pytest.mark.parametrize("module", ["qpdstat", "smartpca"])
def test_eigenstrat_tools_require_explicit_sheep_autosomal_count(tmp_path, chromosome_line, module):
    cfg = parameter_config(tmp_path, module)
    path = Path(cfg["parameter_file"])
    path.write_text(path.read_text().replace("numchrom: 26\n", chromosome_line), encoding="utf-8")
    with pytest.raises(ValueError, match="numchrom: 26"):
        build_plan(module, cfg, tmp_path / "run")


def test_execute_argv_paths_with_spaces_and_metacharacters(tmp_path):
    # A synthetic test double, not a PSMC binary. Any shell interpolation would
    # split the input path and/or create an unwanted marker file.
    program = Path(write_file(tmp_path / "fake psmc", "#!/usr/bin/env python3\nimport json, pathlib, sys\npathlib.Path('result.psmc').write_text(json.dumps(sys.argv[1:]))\n"))
    program.chmod(0o755)
    data = write_file(tmp_path / "sample ; touch PWNED ; $(id).psmcfa")
    directory = tmp_path / "run with spaces"
    plan = build_plan("psmc", {"psmcfa": data, "executable": str(program)}, directory)
    result = execute_plan(plan)
    assert result["status"] == "completed"
    assert json.loads((directory / "result.psmc").read_text())[-1] == data
    assert not (directory / "PWNED").exists()
    assert not (tmp_path / "PWNED").exists()
    assert (directory / "run.json").is_file()
    with pytest.raises(FileExistsError):
        execute_plan(plan)


def test_execute_explicitly_disables_shell_and_detects_missing_outputs(tmp_path, monkeypatch):
    data = write_file(tmp_path / "sample.psmcfa")
    plan = build_plan("psmc", {"psmcfa": data, "executable": sys.executable}, tmp_path / "run")
    calls = []

    def fake_run(argv, **kwargs):
        calls.append((argv, kwargs))
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(subprocess, "run", fake_run)
    record = execute_plan(plan)
    assert calls[0][1]["shell"] is False
    assert record["status"] == "failed"
    assert record["missing_expected_outputs"] == ["result.psmc"]


def test_existing_directory_and_missing_executable_are_rejected(tmp_path):
    data = write_file(tmp_path / "sample.psmcfa")
    plan = build_plan("psmc", {"psmcfa": data}, tmp_path)
    with pytest.raises(ValueError, match="already exists"):
        execute_plan(plan)
    plan = build_plan("psmc", {"psmcfa": data, "executable": "not-a-real-program-4a7e9d"}, tmp_path / "run")
    with pytest.raises(ValueError, match="executable not found"):
        execute_plan(plan)
    assert not (tmp_path / "run").exists()


def test_default_cli_prints_json_and_does_not_execute(tmp_path, capsys):
    config = write_file(tmp_path / "config.json", json.dumps({"psmc": {"psmcfa": "missing.psmcfa"}}))
    output_dir = tmp_path / "run"
    assert main(["--config", config, "--module", "psmc", "--run-dir", str(output_dir)]) == 0
    printed = json.loads(capsys.readouterr().out)
    assert printed["status"] == "plan_only"
    assert printed["blockers"]
    assert not output_dir.exists()


def test_bundled_examples_plan_all_modules_without_creating_outputs(tmp_path):
    root = Path(__file__).resolve().parents[1]
    config = json.loads((root / "config/legacy.example.json").read_text(encoding="utf-8"))
    for module, settings in config.items():
        plan = build_plan(module, settings, tmp_path / module, root)
        assert plan["status"] == "plan_only"
        assert not (tmp_path / module).exists()
