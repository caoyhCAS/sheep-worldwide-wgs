# Window-statistics reconstruction

`window_stats.py` uses only the Python standard library. It is newly written from
the methods in DOI **10.1093/molbev/msab353**, not recovered author code. Its
validation and explicitly listed statistical choices are reconstruction choices.
This utility summarizes statistics generated elsewhere; it does not calculate
FST, nucleotide diversity, D, or fd from genotypes.

## Normalized input contract

Every window table is a tab-delimited file with a header and unique
`chrom,start,end` keys. Coordinates are **1-based inclusive**: `1 <= start <= end`,
matching VCFtools `BIN_START`/`BIN_END` after header normalization.
Convert native tool coordinates when constructing these normalized tables; this
script never converts coordinates or joins tables by row position. Overlapping
windows are permitted. Duplicate keys, malformed rows, missing columns, infinity,
and nonnumeric data fail. Every original column is preserved; output-column
collisions fail. Output directories are created. Input and output cannot be the
same path.

Explicit missing tokens are `NA`, `NaN`, `-NaN`, `+NaN` (case insensitive), `.`, and empty fields.
They are allowed only for `fd`, or `outliers --missing exclude`; a missing statistic
is not replaced by zero. Empty input tables fail. Output missing values are `NA`.

## PBS

Input columns: `chrom,start,end,fst_cha_efr,fst_cha_hus,fst_efr_hus`.
All three FST estimates **must already be joined on the exact same interval**;
different interval grids must not be joined by row order or approximate overlap.
The single input table prevents independent row-order joins inside this script.
The caller remains responsible for truthful keys when normalizing upstream data.

The focal Changthangi statistic is
`(-log1p(-FST_CHA_EFR) - log1p(-FST_CHA_HUS) + log1p(-FST_EFR_HUS)) / 2`.
FST >= 1 is rejected, not capped. Negative FST estimates fail by default;
`--negative-fst zero` explicitly truncates them to zero before transformation.
The policy is recorded in output. A negative calculated PBS is retained.

```sh
python scripts/window_stats.py pbs --input pairwise.tsv --output results/pbs.tsv
```

## Reduction of diversity (ROD)

Input columns: `chrom,start,end,pi_landrace,pi_improved`.
ROD is `1 - pi_improved/pi_landrace`. Both inputs must be finite, nonnegative,
and calculated using matching sites/windows and units. A zero landrace denominator
produces `rod=NA, rod_status=zero_denominator`. Negative ROD is valid and retained.

```sh
python scripts/window_stats.py rod --input diversity.tsv --output results/rod.tsv
```

## fd cleaning, Z scores, and BH

Input columns: `chrom,start,end,fd,D` (D is uppercase). Finite D must be in [-1,1];
out-of-range D is treated as malformed input and rejected, not silently cleaned.
Finite fd values become zero if fd < 0, D < 0, or fd > 1. The fd > 1 rule is applied
even when D == 0 as a conservative reconstruction decision. If either input is
missing, the cleaned fd, Z, P, and Q remain missing and the window is excluded
from the test family. The original input values and cleaning status are retained.

The paper reports a Z transform and Benjamini–Hochberg correction, but not its
exact centering, standard-deviation convention, tail, or family definition.
Therefore this reconstruction **requires** `--tail upper --sd ddof1` to accept:

- Center over all nonmissing cleaned windows in this one input table.
- Scale by their sample standard deviation (`n-1` denominator).
- Use upper-tail standard-normal P = `erfc(Z/sqrt(2))/2`.
- Apply BH once over precisely that set of nonmissing P values.

Put one scientifically intended test family (for example one recipient/donor
comparison) in each input table. Do not silently pool unrelated comparisons.
With fewer than two usable windows or zero variance, all test statistics remain
NA with an explanatory status. This Z calibration is an exploratory reconstruction,
**not a validated null model**. fd is bounded and neighboring windows overlap;
normality and BH dependence assumptions are not established by this script. A
small adjusted P is not by itself proof of introgression.

```sh
python scripts/window_stats.py fd --input fd.tsv --tail upper --sd ddof1 --output results/fd.tsv
```

The paper conflicts on fd window steps: prose says 50 kb, whereas the command
`-w 100000 -m 500 -s 20000` says 20 kb. This helper does not choose a step or
recalculate windows; record the upstream choice separately.

## Empirical high-tail outliers

The paper uses top 1% FST/XP-CLR and 99th-percentile PBS, and top 5% CLR. The exact
quantile algorithm and boundary-tie rule were not given. This implementation uses
linear-interpolated **type-7** sample quantiles, at `1 - top_fraction`, and selects
all values **>=** the threshold. Boundary ties can yield more than the nominal
fraction, even all windows when every score is equal. This is not a randomization
P value and does not perform the paper's separate SV permutation procedure.
Missing input fails unless `--missing exclude`; all-missing input always fails.
The threshold and finite window count are recorded in each output row.

```sh
python scripts/window_stats.py outliers --input pbs.tsv --column pbs --top-fraction 0.01 --output results/pbs_outliers.tsv
```

## Incomplete lineage sorting (ILS)

For total branch length `t` in generations, mean ancestral tract length
`L = 1/(r*t)` and `x = m*r*t`, gamma shape-2 survival is
`P(length > m) = exp(-x)*(1+x)`. We compute `log(P) = log1p(x)-x` to retain a useful
result even when P underflows. `probability_status=underflow` means the reported
zero is a floating-point limitation, **not an exactly zero ILS probability**.

Here `t = branch_factor * divergence_years / generation_years`; the default
`branch_factor=2` assumes the sum of two equally long descendant branches. The
paper describes branch length but does not resolve every conversion convention.
Check the divergence event, total branch definition, and local recombination
model before interpreting results. Rate and time are required, not inferred.
This simple tract model does not account for ancestral recombination or population
structure beyond the stated gamma approximation.

The following is a **demonstration parameter set**, not an assertion that these
are the paper's HBB inputs:

```sh
python scripts/window_stats.py ils --tract-bp 85400 --divergence-years 11000 --generation-years 3 --recombination-rate 1.5e-8 --branch-factor 2 --output results/demo/ils.tsv
```
