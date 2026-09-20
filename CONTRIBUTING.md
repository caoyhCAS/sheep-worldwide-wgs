# Contributing

Keep original-source facts separate from reconstruction choices. Every new
analysis needs a paper/method reference, explicit input/output schema,
parameter provenance, and tests. Do not substitute an undocumented default
for a missing paper parameter without labeling it and requiring review.

Run `python -m pytest -q`, the Snakemake synthetic demo and
`python scripts/check_demo.py` before proposing changes. Never commit real
large sequencing files, credentials or third-party code without permission.

A passing synthetic CI is not validation of a scientific conclusion. Changes
to estimators, sample sets, masks, coordinate systems or test families need
scientific review as well as software tests.
