"""Small deterministic checks for the reconstructed statistics (no genomics binaries)."""

import csv
import importlib.util
import math
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "window_stats.py"
SPEC = importlib.util.spec_from_file_location("window_stats", SCRIPT)
stats = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stats)


def table(tmp_path, text):
    source = tmp_path / "input.tsv"
    source.write_text(text, encoding="utf-8")
    return source


def test_pbs_known_value():
    assert stats.pbs(0.2, 0.3, 0.1) == pytest.approx((-math.log(0.8) - math.log(0.7) + math.log(0.9)) / 2)


def test_pbs_preserves_negative_final_value():
    assert stats.pbs(0, 0, 0.5) < 0


def test_pbs_negative_default_error():
    with pytest.raises(ValueError, match="negative FST"):
        stats.pbs(-0.1, 0.2, 0.3)


def test_pbs_explicit_negative_truncation():
    assert stats.pbs(-0.1, 0.2, 0.3, negative_fst="zero") == stats.pbs(0, 0.2, 0.3)


@pytest.mark.parametrize("invalid", [1, 1.1, float("inf"), float("nan")])
def test_pbs_invalid_fst(invalid):
    with pytest.raises(ValueError):
        stats.pbs(invalid, 0.2, 0.3)


def test_pbs_small_fst_precision():
    assert stats.pbs(1e-20, 1e-20, 0) == pytest.approx(1e-20, rel=1e-12, abs=0)


def test_rod_zero_denominator():
    assert stats.rod(0, 0.1) is None
    assert stats.rod(0, 0) is None


def test_rod_values_and_negative_result():
    assert stats.rod(0.1, 0.05) == pytest.approx(0.5)
    assert stats.rod(0.1, 0.2) == pytest.approx(-1)


@pytest.mark.parametrize("pair", [(-0.1, 0.2), (0.1, float("nan")), (float("inf"), 0.1)])
def test_rod_invalid_values(pair):
    with pytest.raises(ValueError):
        stats.rod(*pair)


@pytest.mark.parametrize("fd,d", [(-0.1, 0.1), (0.2, -0.1), (1.1, 0.1), (1.1, 0)])
def test_fd_cleanup(fd, d):
    assert stats.clean_fd(fd, d) == (0, "set_to_zero")


def test_fd_missing_not_zero():
    assert stats.clean_fd(None, 0.2) == (None, "missing")
    assert stats.clean_fd(0.2, None) == (None, "missing")
    assert stats.clean_fd(1, 0.1) == (1, "retained")


@pytest.mark.parametrize("d", [-1.01, 1.01, float("inf"), float("nan")])
def test_fd_rejects_malformed_d(d):
    with pytest.raises(ValueError):
        stats.clean_fd(0.1, d)
    with pytest.raises(ValueError):
        stats.clean_fd(None, d)


def test_fd_z_and_bh_missing_excluded():
    z, p, q, status = stats.fd_z_bh([0, 0.5, None, 1])
    assert status == "ok"
    assert z == [-1, 0, None, 1]
    assert p[1] == 0.5
    assert q[2] is None
    assert q[3] == pytest.approx(3 * 0.5 * math.erfc(1 / math.sqrt(2)))
    assert q[1] == pytest.approx(0.75)


@pytest.mark.parametrize("values,status", [([None], "insufficient_windows"), ([0.1, None], "insufficient_windows"), ([0, 0], "zero_variance")])
def test_fd_undefined_variance(values, status):
    z, p, q, actual = stats.fd_z_bh(values)
    assert actual == status
    assert z == p == q == [None] * len(values)


def test_bh_known_example_and_ties():
    result = stats.benjamini_hochberg([0.01, 0.04, None, 0.03, 0.002])
    assert result[:2] == pytest.approx([0.02, 0.04])
    assert result[2] is None
    assert result[3:] == pytest.approx([0.04, 0.008])
    assert stats.benjamini_hochberg([0.01, 0.01, 1]) == pytest.approx([0.015, 0.015, 1])


def test_bh_invalid_probability():
    with pytest.raises(ValueError):
        stats.benjamini_hochberg([1.01])


def test_quantile_linear_type7():
    assert stats.empirical_threshold([0, 1, 2, 3, 4], 0.25) == 3
    assert stats.empirical_threshold([0, 10], 0.01) == pytest.approx(9.9)


def test_quantile_ties_include_all():
    values = [1, 1, 1, 1]
    threshold = stats.empirical_threshold(values, 0.01)
    assert sum(value >= threshold for value in values) == 4


def test_quantile_nonfinite_and_empty_fail():
    for values in ([], [1, float("nan")], [float("inf")]):
        with pytest.raises(ValueError):
            stats.empirical_threshold(values, 0.01)


def test_ils_analytic_formula_and_branch_factor():
    result = stats.ils_probability(10000, 10000, 3, 1e-8)
    x = 10000 * 1e-8 * (2 * 10000 / 3)
    assert result["probability"] == pytest.approx((1 + x) * math.exp(-x))
    assert result["branch_generations"] == pytest.approx(2 * 10000 / 3)
    assert stats.ils_probability(0, 10000, 3, 1e-8)["probability"] == 1
    assert stats.ils_probability(10000, 10000, 3, 1e-8, 1)["probability"] > result["probability"]


def test_ils_underflow_keeps_finite_log():
    result = stats.ils_probability(1e8, 15000, 3, 1e-8)
    assert result["probability"] == 0
    assert result["probability_status"] == "underflow"
    assert math.isfinite(result["log_probability"])
    assert result["log_probability"] < -9000


@pytest.mark.parametrize("values", [(-1, 10, 3, 1e-8), (1, 0, 3, 1e-8), (1, 10, 0, 1e-8), (1, 10, 3, -1e-8)])
def test_ils_invalid_arguments(values):
    with pytest.raises(ValueError):
        stats.ils_probability(*values)


@pytest.mark.parametrize("coordinates", ["-1\t10", "0\t10", "11\t10", "0.5\t10"])
def test_bad_intervals(tmp_path, coordinates):
    source = table(tmp_path, f"chrom\tstart\tend\tscore\n1\t{coordinates}\t1\n")
    with pytest.raises(ValueError):
        stats.read_windows(source, ["score"])


def test_duplicate_keys_rejected(tmp_path):
    source = table(tmp_path, "chrom\tstart\tend\tscore\n1\t1\t10\t1\n1\t01\t010\t2\n")
    with pytest.raises(ValueError, match="duplicate window"):
        stats.read_windows(source, ["score"])


def test_single_site_inclusive_window_is_valid(tmp_path):
    source = table(tmp_path, "chrom\tstart\tend\tscore\n1\t1\t1\t0.5\n")
    _, rows = stats.read_windows(source, ["score"])
    assert len(rows) == 1


def test_missing_column_and_malformed_rows(tmp_path):
    source = table(tmp_path, "chrom\tstart\tend\n1\t1\t10\n")
    with pytest.raises(ValueError, match="missing columns"):
        stats.read_windows(source, ["score"])
    source = table(tmp_path, "chrom\tstart\tend\tscore\n1\t1\t10\n")
    with pytest.raises(ValueError, match="wrong number"):
        stats.read_windows(source, ["score"])


def test_explicit_missing_and_infinity_policy():
    assert stats.number("NaN", "fd", missing=True) is None
    assert stats.number("-nan", "fd", missing=True) is None
    assert stats.number("+NaN", "fd", missing=True) is None
    for value in ("NaN", "inf", "hello"):
        with pytest.raises(ValueError):
            stats.number(value, "score")
    with pytest.raises(ValueError):
        stats.number("inf", "score", missing=True)


def test_cli_fd_requires_explicit_statistical_choices(tmp_path):
    completed = subprocess.run([sys.executable, str(SCRIPT), "fd", "--input", "x", "--output", str(tmp_path / "output")], capture_output=True, text=True)
    assert completed.returncode != 0
    assert "--tail" in completed.stderr


def test_cli_fd_preserves_missing_and_writes_parent(tmp_path):
    source = table(tmp_path, "chrom\tstart\tend\tfd\tD\n1\t1\t10\t0\t1\n1\t11\t20\tNA\t1\n1\t21\t30\t1\t1\n")
    output = tmp_path / "nested" / "out.tsv"
    assert stats.main(["fd", "--input", str(source), "--output", str(output), "--tail", "upper", "--sd", "ddof1"]) == 0
    with output.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[1]["fd_clean"] == rows[1]["q_bh"] == "NA"
    assert {row["test_n"] for row in rows} == {"2"}


def test_cli_missing_outliers_are_opt_in(tmp_path):
    source = table(tmp_path, "chrom\tstart\tend\tscore\n1\t1\t10\tNA\n1\t11\t20\t1\n")
    output = tmp_path / "out.tsv"
    args = ["outliers", "--input", str(source), "--output", str(output), "--column", "score", "--top-fraction", "0.01"]
    with pytest.raises(SystemExit):
        stats.main(args)
    stats.main(args + ["--missing", "exclude"])
    with output.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["is_outlier"] == "NA"
    assert rows[1]["is_outlier"] == "1"


def test_cli_rejects_inplace_output(tmp_path):
    source = table(tmp_path, "chrom\tstart\tend\tpi_landrace\tpi_improved\n1\t1\t10\t1\t0.5\n")
    with pytest.raises(SystemExit):
        stats.main(["rod", "--input", str(source), "--output", str(source)])


def test_output_column_collision():
    with pytest.raises(ValueError, match="already exist"):
        stats.append_fields(["pbs"], ["pbs"])
