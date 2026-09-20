# Validation record for reconstruction v0.1.0

Validated locally on 2026-09-20 using Linux, Python 3.12.14,
Snakemake 8.25.5 and pytest 8.3.5.

| Check | Result |
| --- | --- |
| Python unit/integration tests | 99 passed, none skipped |
| Python compilation | Passed for scripts, workflow helpers and tests |
| All real workflow rule types | 28 rule types formatted in a synthetic-input dry-run, including optional ROH, LD and global FST |
| Synthetic demonstration | All 7 jobs completed, including final aggregation |
| Numerical demo assertions | Passed for PBS outliers, ROD missing denominator, fd cleanup/BH missingness and ILS probability |
| Citation metadata | Valid CFF 1.2.0; preferred paper citation contains all 56 listed authors |
| Config JSON/YAML and local documentation links | Parsed / checked |

Reproduce the main checks from the repository root:

```bash
python -m pytest -q -rs
python -m compileall -q scripts workflow/scripts tests
snakemake --snakefile workflow/Snakefile --configfile config/config.test.yaml --dry-run --cores 1
snakemake --snakefile workflow/Snakefile --configfile config/config.test.yaml --cores 1
python scripts/check_demo.py
```

The real-module DAG test creates temporary placeholder paths. It does not
pretend they contain real sequencing data: it checks dependency resolution and
shell command construction, **not** external program execution. Legacy-adapter
execution tests use synthetic Python test doubles, not installed historical
bioinformatics tools.

Not performed: full Conda environment solves, FASTQ/VCF analysis on real samples,
running every external binary, original-data numerical comparison, complete
reproduction of the paper, remote GitHub Actions execution. The included CI
configuration is ready for a repository, but local tests are not reported as
a successful remote CI run.
