# Data and input contracts

## Study resources

Reference: **Oar_rambouillet_v1.0**, NCBI assembly **GCF_002742125.1**.
Do not substitute Oar_v3.1 or a newer assembly while keeping the paper's
coordinates. Paper IRF2BP2 interval: chr25:7,067,974–7,071,785 (1-based).
FASTA, annotations, indexes, contig names and variant coordinates must agree.

The paper examined 810 sheep (72 wild, 738 domestic). Its data-availability
statement says newly generated WGS can be requested from authors for research.
Public sample accessions are in Table S1; RNA and miRNA lists are in S45/S47;
demographic-model sample choices are in S50; wool groups are in S32.
[Download the original supplement](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8826587/supplementaryFiles).
Missing or malformed accession entries are unresolved, not synthesized here.

Do not upload raw datasets or personal credentials into this code repository.
The ignored local `data/`, `resources/`, `external/`, `results/` and `logs/`
directories are working locations, not a substitute for research data storage.

## WGS and population workflow inputs

Use the exact schemas in [workflow/README.md](../workflow/README.md).
Sample IDs are unique. Samples at depth exactly 10 do not satisfy the
paper's >10 threshold. The WGS cohort is explicitly staged using the
manifest depth, then checked against depth measured by this workflow.
Coverage denominators depend on reference contigs/masks; record those choices.

Population inputs require a PASS biallelic SNP VCF and a deliberate contig map
to sheep autosomes 1–26. Excluded or unrecognised contigs must be investigated.
Unrelated IDs must be curated after kinship analysis; this repository does not
pretend that arbitrary tie-breaking reconstructs the exact 714 retained sheep.

Separate unrelated population lists, wool phenotype groups, ancestry-model
cohorts and high-depth demographic groups. They are not interchangeable.
Do not perform SFS/PSMC/SMC++ from a population-structure MAF-pruned VCF.

## Normalized windows

Internal statistics TSVs use unique `chrom`, `start`, `end` keys with
**1-based inclusive** coordinates. Do not provide BED start coordinates.
Headers are case-sensitive. Numeric chromosomes need not be prefixed `chr`.

`prepare_windows.py` joins native VCFtools tab-separated files by
`CHROM,BIN_START,BIN_END`, not row order. Default `--join exact` rejects
nonidentical window sets. `--join intersection` must be explicitly requested;
the script reports input/output counts. Missing numeric cells become `NA`.
The PBS/ROD statistics helpers deliberately reject missing inputs rather than
silently dropping windows; resolve the missingness policy before running them.

For PBS, choose VCFtools `--fst-column MEAN_FST` or `WEIGHTED_FST` explicitly.
This adapter's FST estimator is **not** claimed to be the original paper's R
estimator. Exact estimator and samples must be consistent across all pairs.

The fd helper expects normalized TSV, not the native comma-separated
genomics_general output. Map native `scaffold,start,end,D,fd` to
`chrom,start,end,D,fd`, preserving the numeric values and window coordinates.
Apply the Z/BH helper to one prespecified test family (normally a fixed donor,
recipient and outgroup comparison); concatenating unrelated comparisons changes
the inference. With overlapping windows, independence cannot be assumed.

VCFtools π on SNP-only VCF treats the denominator differently from an analysis
with full callable-site masks. Record missingness, callable regions and variant
ascertainment; a raw window output alone is not proof of an unbiased estimate.

## Deterministic sample helper

`select_samples.py` accepts `sample,population,mean_depth` as tab-separated
columns. It selects three samples per population by default:

```bash
python scripts/select_samples.py --input samples.tsv --mode highest-depth --output psmc_candidates.tsv
python scripts/select_samples.py --input samples.tsv --mode random --seed 2026 --output ld_candidates.tsv
```

The seed is a new recorded reconstruction choice, not the study seed.
Highest-depth ties are broken lexically by sample ID. Too-small groups fail.
The output is a candidate metadata table, **not** a VCFtools keep file or a
PLINK two-column keep file; convert to the required schema deliberately.

## Provenance

Save edited configuration, external tool versions, reference/checksum,
command plans and logs. `record_provenance.py` hashes only named files and
does not inspect environment variables or credentials:

```bash
python scripts/record_provenance.py --input config/config.local.yaml \
  --input samples.tsv --output results/run.inputs.json
```

For large inputs, hashing reads the complete file and has I/O cost. Tool
versions must be recorded separately from this input-file manifest.
