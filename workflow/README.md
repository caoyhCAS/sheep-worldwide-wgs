# Executable core: method-based reconstruction

This is newly written workflow code for DOI `10.1093/molbev/msab353`, **not the
historical author scripts**. The DAG is modular. The synthetic `demo` does not
call any genomic tools and cannot validate biological reproduction.

Run all commands from the repository root:

```bash
snakemake --snakefile workflow/Snakefile --configfile config/config.test.yaml --cores 1
snakemake --snakefile workflow/Snakefile --configfile config/config.local.yaml --dry-run --cores 8
snakemake --snakefile workflow/Snakefile --configfile config/config.local.yaml --use-conda --cores 8
```

Copy `config/config.example.yaml` to `config/config.local.yaml`, edit its paths
and unresolved choices, and select `modules: [wgs]`, `[population]`, `[selection]`
or a combination. The example intentionally contains blocking `null` values.
When chaining WGS to other modules, explicitly set `analysis.vcf` to
`results/wgs/cohort.pass.snps.vcf.gz`. Otherwise an existing indexed PASS SNP VCF
is the analysis entry point. Relative paths are relative to the repository root.

## Scope and data contracts

| Module | Executed calculations | Main output |
| --- | --- | --- |
| `demo` | Synthetic statistics and provenance only | `results/demo/provenance.json` |
| `wgs` | Raw FastQC, Trimmomatic, BWA MEM, sorting, duplicate removal, measured depth gate, GVCF calls, joint calling, hard filtering | `results/wgs/cohort.pass.snps.vcf.gz` and `cohort.pass.indels.vcf.gz` |
| `population` | Curated unrelated cohort, PLINK QC/pruning, IBS distance, diagnostic heterozygosity, individual and pooled π; optional resolved ROH and three-individual LD decay | `results/population/` |
| `selection` | Pairwise and optional all-population windowed FST, population π for downstream ROD | `results/selection/` |

Other analyses are covered by standalone scripts, explicit external-tool
handoffs, or documented gaps; they are not implied to be DAG targets here.

### Reference and WGS

Supply the Oar_rambouillet_v1.0 reference FASTA (assembly GCF_002742125.1), its
`.fai`, BWA `.amb/.ann/.bwt/.pac/.sa` indexes, and a GATK sequence dictionary at
`reference.dict`. Prepare and verify these before running; this workflow does
not mutate the original reference to construct indexes. BAM/reference contig
lengths and names must match. Record the original assembly checksum.

`wgs.samples` is a tab-separated table with header:

```text
sample  r1  r2  mean_depth
```

`sample` is unique, starts with a letter/digit, and otherwise contains only
letters, digits, dot, underscore or hyphen. `r1` and `r2` are paired FASTQs.
The example names and depths are placeholders, not real samples. This first
version supports one paired library per individual. Do not concatenate
libraries with conflicting read-group metadata without a reviewed strategy.

`mean_depth` is a finite nonnegative value from a prior coverage assessment.
It declares the staged short-variant cohort: only manifest values **strictly
greater than 10** enter joint calling. All listed samples are aligned and their
depth measured. Each admitted sample must independently pass a measured >10
gate before HaplotypeCaller can run. If it fails, revise the cohort explicitly;
the workflow does not silently change joint-call sample membership mid-run.
Low-depth manifest samples never enter the joint call even if a new measurement
is higher until the manifest is deliberately updated.

Coverage is reconstructed as the arithmetic mean of `samtools depth -aa` over
all reference positions, including zero-coverage bases; the position count must
equal the summed `.fai` lengths. Tool-default base/read filters apply, overlapping
mates are not collapsed, and no additional depth cap is applied by the portable
SAMtools version. **The paper does not report the exact denominator or depth
options.** This choice affects cohort admission and is not an exact historical
coverage implementation.

Explicitly provide `wgs.trimmomatic_steps` as a list of complete operation
tokens. Adaptor identity, quality trimming and read-length thresholds are not
reported; the repository does not fabricate them. Optional adaptor FASTAs named
inside operation tokens must exist and be recorded in run provenance. Set
`wgs.missing_annotation_policy` to `pass` or `fail` after reviewing GATK missing
annotation handling. This is an unresolved choice, not a paper-derived value.

GATK hard-filter expressions follow the reported thresholds. Each threshold is
emitted as its own named filter; rejecting any such flag reconstructs the logical
OR in the paper, but avoids a compound JEXL expression losing a failing threshold
when a different annotation is missing. This is an explicit implementation
choice. Missing-annotation policy applies independently to each expression;
review filter counts and warnings. Duplicates are removed (`REMOVE_DUPLICATES`
true). BQSR is **not added** because no known-sites resource/explicit BQSR method
is given. GATK defaults, not invented sample ploidy or extra read filters, are
used. PASS SNPs are restricted to biallelic sites; indels are not additionally
restricted to biallelic sites. Whole-cohort `CombineGVCFs` is the paper-described
method and may be expensive at 810-sample scale; no benchmark or sharding claim
is made. No SV workflow is wired to these BAMs automatically.

### Analysis VCF, populations and unrelated samples

`analysis.vcf` is a BGZF VCF with `.tbi` index. Only records with literal `PASS`
filter status and exactly two alleles that are SNPs enter downstream analysis.
Unfiltered `.` status is not treated as PASS. `analysis.contig_map` contains
two tab-separated columns, old and new contig labels, without a header.
Use identity mappings if needed. The explicit `analysis.autosomes` list selects
numeric sheep autosomes in `1..26` after mapping. The pipeline never assumes
RefSeq accession numbers are numeric chromosomes. Check all intended contigs
are retained and do not merge different original contigs under one new label.
SNP IDs are set to `CHROM:POS:REF:ALT` for PLINK bookkeeping. Multi-record
duplicate alleles should be resolved before use; this code does not silently
deduplicate them.

`analysis.populations` is a TSV with `population` and `keep` columns, optionally
`ld_keep`. Each `keep` file has one VCF sample ID per line, no header. Populations
must be nonempty, contain no duplicate IDs and be mutually disjoint. The
validation step rejects unknown members. For `population` runs, all population
members must also be in `population.unrelated_keep`.

`population.unrelated_keep` contains two whitespace-separated columns without a
header: `FID IID`, where both fields equal the VCF sample ID (`--double-id`).
Supply a curated cohort after reviewing KING-robust kinship >0.0884. The paper
does not specify a removal tie-breaker, so this DAG does not choose which
relative to discard or substitute a different kinship estimator.

PLINK MAF ≥0.05, missingness ≤0.10, HWE p ≥0.001 and pruning `50 5 0.2` apply
to the population-structure branch. `--distance square 1-ibs` yields the
proportion-distance `.mdist` matrix. The square shape is a reconstruction output
choice; the paper does not specify matrix shape. HWE is applied jointly to the
provided unrelated cohort, an implementation interpretation requiring review
for mixed populations/species. The paper's reported retained counts are not
hardcoded as expected outputs.

The diversity branch uses the unrelated high-quality calls **before** MAF/HWE
filtering and LD pruning. Individual π is computed in nonoverlapping 1 Mb
VCFtools windows for every individual in the curated unrelated list, in
`results/population/individual_pi/`. The same computation pooled within each
manifest population is a supplementary output in `results/population/pi/`, not
a substitute for the individual statistic. The paper's windowed expected
heterozygosity (`He`) implementation remains unresolved. The supplementary PLINK
`.het` report is a genome-wide diagnostic homozygosity/inbreeding report, **not
the paper's windowed He estimate**. Callability/missing-site treatment can bias
π from a SNP-only VCF and must be assessed; no callable-base mask is inferred.

ROH is disabled until explicitly resolved: the paper lists
`--homozyg-window-kb 200`, which is not documented in the current
[PLINK 1.9 ROH interface](https://www.cog-genomics.org/plink/1.9/ibd).
`--homozyg-kb 200` instead changes minimum segment length and is not silently
substituted. Enabling ROH requires a written `resolution_note` and a nonempty
`resolved_args` list. The unambiguous density 50, window het 1 and window SNP 50
arguments are already supplied; the resolved list adds the user's decision.
Other PLINK defaults remain version dependent. No FROH denominator is invented.

LD decay is disabled by default. Enabling it requires a `ld_keep` file containing
exactly three members per population. Prepare these using a recorded random
seed/selection rule. PopLDdecay `-MaxDist 300` is in kb (300 kb), while its other
filters remain tool defaults. It is applied to the unpruned PASS SNP VCF.

### Selection

Comparisons map unique names to two population names. The analysis VCF is not
subjected to the population-structure MAF/HWE/LD-pruning filters for selection.
The paper reports 50 kb windows in both sections but a 10 kb step in Methods
versus 25 kb step in Results. `fst_step_bp` is mandatory; enter the chosen step
and a `fst_step_resolution` explanation. `fst_window_bp` is explicit (the example
uses 50 kb). With `global_fst: true`, VCFtools also computes an all-population
scan using every manifest population, corresponding to the reported scan across
158 populations. A different population set is a deliberate reconstruction
choice and must be recorded.
VCFtools output contains both weighted and mean FST; use the **mean FST** column
for a reconstruction of the reported average-FST selection ranking, and record
your treatment of undefined/negative estimates and ties.

Population π in `results/selection/pi/` uses nonoverlapping 10 kb windows. Join
these by chromosome/start/end before downstream `1 - π_improved/π_landrace`
calculation. The selection module produces statistics, not a completed candidate
gene list; use the documented scripts to calculate ROD and empirical outliers.

## Software and validation limits

`environment.yml` pins the workflow runner. Rule-specific `envs/*.yaml` are
**portable environment specifications**, not fully solved lockfiles and not an
exact historical environment. Older paper versions can differ numerically.
The test suite does not claim that all Conda specifications solve on every
platform, or that all bioinformatics binaries have been executed on real WGS.
Use Linux, inspect Conda solves before scheduling a large run, retain resolved
package metadata and compare small real-data subsets before scaling up.

| Program | Paper | Portable specification |
| --- | --- | --- |
| FastQC | 0.11.9 | 0.11.9 |
| Trimmomatic | 0.36 | 0.39 |
| BWA | 0.7.8 | 0.7.17 |
| SAMtools | 1.9 | 1.20 |
| GATK | 4.1.2.0 | 4.5.0.0 |
| PLINK | 1.90/1.90p | 1.90b6.21 |
| VCFtools | 0.1.17 / 1.17 as written | 0.1.16 |
| PopLDdecay | 3.41 | 3.42 |
| BCFtools | Not reported for this preprocessing | 1.20, explicit portability helper |

Official interface references used to check syntax:
[GATK SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360037434011-SelectVariants),
[SAMtools depth](https://www.htslib.org/doc/samtools-depth.html),
[PLINK distance matrices](https://www.cog-genomics.org/plink/1.9/distance), and
[PopLDdecay](https://github.com/hewm2008/PopLDdecay).
