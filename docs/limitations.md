# Reconstruction limits and unresolved choices

**方法重建版，非原始代码 / Method-based reconstruction, not original code.**

Primary evidence: [paper](https://doi.org/10.1093/molbev/msab353),
[full text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8826587/fullTextXML),
[supplementary archive](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8826587/supplementaryFiles).
The examined supplement contains figure information and Tables S1–S50, not
the authors' complete pipeline. Absence from these examined materials is not
a claim that no original scripts exist elsewhere.

## Decisions that must not be hidden

| Issue | Evidence / unresolved detail | Reconstruction treatment |
| --- | --- | --- |
| fd window step | Prose: 50 kb; explicit command: `-s 20000` | Require 20 kb or 50 kb explicitly; preserve choice in command plan |
| Global selection FST step | Methods: 50 kb windows / 10 kb step; Results: 25 kb step | Explicit setting and resolution note; neither is silently asserted as historical truth |
| ROH flag | Paper reports `--homozyg-window-kb 200`, absent in current PLINK 1.9 manual | Disabled until user resolves flag/version; `--homozyg-kb` is not an automatic equivalent |
| ILS branch length | Reported expected lengths imply two branches, `t=2*years/generation` | Expose branch factor, preserve log probability, document discrepancy below |
| SNP sets | 24,363,608 post-filter/pruning SNPs vs a distinct 220,350 unlinked set | Do not equate these datasets or promise identical counts |
| Trimming | Adapter sequences, quality settings not supplied | User-provided Trimmomatic arguments required |
| Coverage | Whole-genome depth definition and masks not fully supplied | Explicit cohort and measured depth gate; this implementation's denominator is documented |
| SV consensus | Caller matching tolerances and reconciliation script unavailable | Method specification only; no fabricated consensus algorithm |
| fastsimcoal | `.tpl`, `.est`, ancestral calls, masks and prior ranges absent | Caller-supplied vetted inputs required; wrapper does not invent models |
| PSMC | Consensus conversion filters and exact sample choices missing | Accept prepared `.psmcfa`; no SNP-only VCF-to-consensus shortcut |
| SMC++ | Callable/SNPable masks absent | Accept prepared masked `.smc.gz`; do not treat absent VCF records as known invariant sites |
| fd significance | Z-transform tail, SD convention, multiple-testing scope not fully specified | New statistical policy explicit in helper CLI/documentation; not validated original inference |
| qpDstat sheep chromosomes | Official ADMIXTOOLS v7.0.1 defaults to 22 autosomes and does not read `numchrom` from the parameter file | Execution gated until a patched or newer build's chromosome-26 support is verified and documented; setting `numchrom: 26` alone is insufficient for stock v7.0.1 |
| sNMF | Seeds/repetitions not reported | Manual handoff; no substitution of ADMIXTURE or newer LEA as the historical executable |
| Gene enrichment | Historical DAVID 6.8 background and database state unavailable | Manual handoff with saved gene universe and results; no current database asserted identical |
| RNA / miRNA | Adapter/quality details, several tool versions and databases missing | Documented handoff; no full reproduction claimed |

## ILS numerical discrepancy

For a gamma survival function with shape 2, `x = tract_length / L`,
`log(P) = log1p(x) - x`, and `L = 1 / (r * t)`.
At `r=1.5e-8` per bp per generation, generation time 3 years,
tract length 85,400 bp and `t=2*divergence_years/3`:

- Asiatic mouflon at 11 ka: probability approximately 0.000865.
- European mouflon at 5–6 ka: approximately 0.0737–0.0365.
- Argali at 2.36 Ma and urial at 1.26 Ma: floating-point survival can
  underflow; a finite log probability must be retained.

These evaluations of the reported formula are not all zero, unlike the
paper's wording. The helper reports the formula's result, not the published
zero assertion. This is a computational discrepancy, not by itself a revised
biological conclusion; divergence and recombination assumptions matter.

## Historical qpDstat version issue

Inspection of official [v7.0.1 globals.h](https://github.com/DReichLab/AdmixTools/blob/v7.0.1/src/globals.h)
and [qpDstat.c](https://github.com/DReichLab/AdmixTools/blob/v7.0.1/src/qpDstat.c)
found a human-default `numchrom=22`, chromosome filtering, and no parameter-file
reader for `numchrom` in that version. The ordinary sheep autosomes 23–26 can
therefore be omitted by an unmodified historical executable. Later source
supports the override. This audit does **not** establish which local changes
the study authors made or which chromosomes their analysis actually used.
The reconstruction must not claim that a parameter override fixes the stock
historical build. A verified patched/newer binary is a recorded version
departure; it is required before execution through this adapter.

## Data and validation limits

The paper states that newly generated WGS are available from authors upon
request for research purposes. Other public samples are listed in Table S1.
The source table has missing and apparently non-accession entries; accession
numbers must be checked before downloading. No data requests are sent and no
large sequence downloads are started by these scripts.

Synthetic tests check mathematics, schemas, command construction, and workflow
wiring. They do **not** test all legacy executables, recreate the 810-sample
study, reproduce its figures, or establish biological validity. The default
demo runs only the repository's Python code, not WGS tools.
