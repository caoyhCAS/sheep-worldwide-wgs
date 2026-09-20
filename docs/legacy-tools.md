# Prepared-input adapters for legacy methods

These are newly written **method-reconstruction adapters, not the authors' original scripts**. They neither download external software nor claim to reproduce the paper's biological results. Seven adapters are available in `scripts/legacy_tools.py`; other methods remain documented handoffs.

## Plan first, execute explicitly

From the repository root, inspect a non-mutating JSON plan:

```bash
python scripts/legacy_tools.py --config config/legacy.example.json --module psmc --run-dir results/legacy/psmc_sample_01
```

The plan lists the exact argument vector, expected outputs, required input files, missing-file blockers, and the reported/requested tool version. Missing inputs are shown without creating directories or running the program. The example paths and populations are placeholders; edit them for your own reviewed inputs. Paths in JSON are relative to the **invocation working directory**, not the JSON file. Input paths in EIGENSOFT/ADMIXTOOLS parameter files are relative to the parameter file's directory.

After reviewing the plan and preparing the inputs, append `--execute` to the same command. Each execution requires a **new, nonexistent run directory**; there is no overwrite/force option. Arguments are passed with `shell=False`; spaces and shell metacharacters in paths are not evaluated by a shell. No arbitrary `extra_args` field is accepted. Only the configured external executable/script runs; install trusted releases yourself.

An execution saves `run.json`, standard output and standard error. The record says explicitly that the requested version has **not** been verified automatically. Record your executable version/build and the external repository commit separately before analysis. Exit code zero plus expected-file existence is a technical check, not statistical validation. Failed-run directories are retained for inspection and cannot be silently reused.

## What each adapter does

| Module | Required prepared input | Program/parameters | Expected scientific output inside run directory |
|---|---|---|---|
| `psmc` | Callable-mask-aware diploid `.psmcfa` | `psmc -N25 -t15 -r5 -p 4+25*2+4+6` | `result.psmc`, in native scaled units |
| `smcpp` | Explicit list of prepared SMC++ files | `smc++ estimate -o estimate 1.51e-8 ...` | `estimate/model.final.json`; not a years-scaled plot |
| `fastsimcoal` | User-reviewed `.tpl`, `.est`, unfolded multidimensional `.obs`, explicit seed | One fsc26 run: `-n 1000000 -d -M -l 25 -L 65 --multiSFS --seed ...` | `model/` results directory; exact contents depend on build |
| `qpdstat` | Restricted parameter file, EIGENSTRAT/quartet files, verified chromosome-support acknowledgement | `qpDstat -p analysis.par`; stock v7.0.1 execution blocked | `results.txt` (program stdout, including D and Z statistics) |
| `abba` | External script checkout, genomics-general genotype file and population assignments | `ABBABABAwindows.py`, 100 kb windows, 500 minimum sites, user-selected step | `windows.csv` with raw D/fd statistics |
| `popgen` | Same input contract; exactly two population names | `popgenWindows.py --analysis popDist popPairDist`, same window settings | `windows.csv` with within/between-population statistics |
| `smartpca` | Restricted parameter file and prepared EIGENSTRAT files | `smartpca -p analysis.par` | `results.evec`, `results.eval` |

PSMC consensus/mask construction and selection of the three highest-depth individuals are upstream tasks. A SNP-only input does not encode callable invariant sites. Use the study's mutation rate, `1.51e-8` per site per generation, and 3-year generation time only at the appropriate scaling/plotting stage; this wrapper does not relabel raw PSMC or SMC++ output as calendar years. CLI references: [PSMC](https://github.com/lh3/psmc/blob/master/README), [SMC++](https://github.com/popgenmethods/smcpp).

For fsc26, the wrapper stages copies as `model.tpl`, `model.est`, and `model_DSFS.obs`. Supply the correctly polarized multidimensional spectrum, callable monomorphic count, population order, demographic events, mutation model/rate, parameter ranges and units yourself. The `.est` must include a `reference` parameter for the reported `-l 25` switch: in v2.6 this is the number of initial cycles retaining monomorphic-site information, not a minimum number of cycles. A keyword check catches omission but does not validate model biology. The seed range is 1–1,000,000. This adapter executes **one** optimization; orchestrate the paper's 100 independent restarts using distinct recorded seeds and new run directories. It does not manufacture the three model specifications, block-bootstrap data or 100 × 20 bootstrap restarts. See the [fsc26 manual](https://cmpg.unibe.ch/software/fastsimcoal26/man/fastsimcoal26.pdf), especially §§7–8.

## Introgression ambiguity and callable sites

`step_bp` is mandatory and accepts only `20000` or `50000`. The paper's prose says 50 kb whereas the printed command says `-s 20000`; the example JSON explicitly chooses the printed-command value, **not** a resolved historical truth. Run sensitivity analyses or obtain author clarification before combining outputs. Both adapters use the same explicit selection.

For `abba`, `populations` is ordered `[P1, P2, P3, O]`: Menz reference, recipient, wild donor, outgroup. The ancestral-state pseudo-outgroup required for some donor comparisons is not synthesized. For `popgen`, supply `[donor, recipient]`, then make a separate configuration/run for `[donor, reference]`. The genotype file must use the chosen `phased`, `pairs`, `haplo` or `diplo` format; it is **not** a VCF. The population file is the upstream two-column sample/population file. Optional `ploidy_file` maps sample IDs to ploidy; its suitability is the caller's responsibility. Absolute diversity/divergence require callable invariant sites and consistent missingness masks, not only segregating variants.

Install the complete external `genomics_general` checkout so `genomics.py` is adjacent to each script. CLI syntax was checked against the upstream [ABBABABAwindows.py](https://github.com/simonhmartin/genomics_general/blob/master/ABBABABAwindows.py) and [popgenWindows.py](https://github.com/simonhmartin/genomics_general/blob/master/popgenWindows.py). The paper's exact revision is unreported. No third-party code is copied into this repository. Raw CSVs are not automatically filtered, converted to BED, or assigned P/FDR values by these adapters.

## Restricted parameter files

For `qpdstat`, required keys are `genotypename`, `snpname`, `indivname`, `popfilename`, and **`numchrom: 26`**. The quartet file supplies explicit W/X/Y/Z tests. Optional allowed keys are `f4mode` (must be `NO`), `printsd`, `blgsize`, and `inbreed`; other settings retain the installed binary's defaults. No significance-direction reversal is performed. The study's criterion is `|Z| > 3`. See [ADMIXTOOLS D-statistics documentation](https://github.com/DReichLab/AdmixTools/blob/master/README.Dstatistics).

**Critical version issue:** stock [ADMIXTOOLS v7.0.1 `qpDstat.c`](https://github.com/DReichLab/AdmixTools/blob/v7.0.1/src/qpDstat.c) filters chromosomes using `numchrom` but does not parse that parameter. Its [global default is 22](https://github.com/DReichLab/AdmixTools/blob/v7.0.1/src/globals.h), so an unmodified binary drops sheep chromosomes 23–26 even when `numchrom: 26` appears in the file. [Current upstream source](https://github.com/DReichLab/AdmixTools/blob/master/src/qpDstat.c) does read the parameter; that does not establish compatibility of an installed build. Accordingly, execution is blocked by default. To proceed, independently audit a patched legacy or appropriate newer binary, set `chromosome_26_support_confirmed` to `true`, identify it in `requested_version` (not just stock `7.0.1`), and record the tested patch/build and chromosome-inclusion evidence in `supported_build_evidence`. This is a caller attestation, **not** automated verification, and the patch/version change is an explicit reconstruction departure. Do not renumber or discard chromosomes silently to evade this check.

For `smartpca`, required keys are `genotypename`, `snpname`, `indivname`, and **`numchrom: 26`**. This explicitly adapts the human-oriented default to sheep autosomal inputs labeled 1–26; the original paper does not report this conversion detail. Output names are normalized to `results.evec` and `results.eval`; source parameter files are unchanged. Supported optional scalar keys are `numoutevec`, `numoutlieriter`, `numoutlierevec`, `outliersigmathresh`, `usenorm`, `altnormstyle`, `missingmode`, and `lsqproject`. Other defaults come from the installed version. See [EIGENSOFT documentation](https://github.com/DReichLab/EIG/blob/master/POPGEN/README).

Only `key: value` lines and full-line `#` comments are accepted. Duplicate keys, unsupported directives, macros, and non-scalar optional values fail closed. Input dependencies are checked, then symlinked under safe local names so original paths with spaces do not enter legacy parameter-file parsers. Extra output directives are rejected to avoid writing outside the new run directory. This deliberately limited format is not a general-purpose interpreter for all ADMIXTOOLS/EIGENSOFT options.

## Explicit handoffs, not automated wrappers

- XP-CLR v1.0: the paper reports `-w1 0.005 200 2000 2 -p0 0.95`. Supply correctly ordered genotype and genetic-map inputs and validate the archived binary's full CLI. The current [hardingnj/xpclr](https://github.com/hardingnj/xpclr) is a different implementation with a different interface; substituting it would not reproduce the reported command. No executable XP-CLR adapter is claimed here.
- sNMF v1.2: K = 1–20 with cross-entropy selection is reported, but input conversion, initialization seeds and repetition count require decisions. No sNMF optimization is launched by these adapters.
- SNeP, SweeD, RNA/small-RNA analysis, SV reconciliation and gene-set enrichment need additional method-specific inputs and review. In particular, a SweeD `-grid` value denotes evaluated grid points, not automatically independent 10 kb windows; no invented command is supplied.
- The three-tool SV genotype reconciliation and fastsimcoal scenario/constraint files are not recoverable simply from the reported software names. They remain separate reconstruction tasks.

The tests validate command construction, required inputs, the window-step conflict, output isolation, and shell-safe argument handling using synthetic test doubles. They do **not** validate legacy program installation, numerical equivalence, full-data execution or biological replication.
