import csv
import pytest

from prepare_windows import join_windows, read_windows, main as prepare
from select_samples import select
from record_provenance import sha256


def test_join_uses_coordinates_not_order():
    a = {("2", 1, 10): ".1", ("1", 1, 10): ".2"}
    b = dict(reversed(list(a.items())))
    assert join_windows([a, b]) == list(a)


def test_join_rejects_missing_windows():
    with pytest.raises(ValueError, match="differ"):
        join_windows([{("1", 1, 10): 1}, {("1", 11, 20): 2}])


def test_intersection_explicit():
    a = {("1", 1, 10): 1, ("1", 11, 20): 2}
    assert join_windows([a, {("1", 11, 20): 4}], "intersection") == [("1", 11, 20)]


@pytest.mark.parametrize("body", ["1\t0\t10\t.1\n", "1\t10\t1\t.1\n", "1\t1\t10\tinf\n", "1\t1\t10\t.1\n1\t1\t10\t.2\n"])
def test_invalid_windows(tmp_path, body):
    path = tmp_path / "pi.tsv"
    path.write_text("CHROM\tBIN_START\tBIN_END\tPI\n" + body)
    with pytest.raises(ValueError):
        read_windows(path, "PI")


def test_native_rod_cli(tmp_path):
    a, b, dest = [tmp_path / p for p in ["a.pi", "b.pi", "out.tsv"]]
    a.write_text("CHROM\tBIN_START\tBIN_END\tPI\n1\t1\t10000\t.2\n")
    b.write_text("CHROM\tBIN_START\tBIN_END\tPI\n1\t1\t10000\t-nan\n")
    prepare(["rod", "--landrace", str(a), "--improved", str(b), "--output", str(dest)])
    with dest.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert rows[0]["pi_improved"] == "NA"
    assert rows[0]["start"] == "1"


def rows(pop="A"):
    return [{"sample": f"{pop}{i}", "population": pop, "mean_depth": str(i+10)} for i in range(5)]


def test_highest_depth():
    assert [r["sample"] for r in select(rows(), "highest-depth", 3)] == ["A4", "A3", "A2"]


def test_random_reproducible_and_order_invariant():
    assert select(rows(), "random", 3, 42) == select(list(reversed(rows())), "random", 3, 42)


def test_random_population_independent():
    chosen = select(rows()+rows("B"), "random", 3, 42)
    assert [r for r in chosen if r["population"] == "A"] == select(rows(), "random", 3, 42)


@pytest.mark.parametrize("data,mode,n,seed", [(rows(), "random", 3, None), (rows(), "highest-depth", 6, None), (rows()+rows()[:1], "highest-depth", 3, None), ([], "highest-depth", 3, None)])
def test_invalid_selection(data, mode, n, seed):
    with pytest.raises(ValueError):
        select(data, mode, n, seed)


def test_hash(tmp_path):
    path = tmp_path / "test.txt"
    path.write_bytes(b"abc")
    assert sha256(path) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
