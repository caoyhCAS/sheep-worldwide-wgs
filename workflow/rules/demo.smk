# Synthetic demonstration: no genomic sequence data or external executables.
rule demo_pbs:
    input:
        table="tests/fixtures/demo_pbs.tsv",
        script="scripts/window_stats.py"
    output:
        "results/demo/pbs.tsv"
    shell:
        "python {input.script:q} pbs --input {input.table:q} --output {output:q}"

rule demo_pbs_outliers:
    input:
        table="results/demo/pbs.tsv",
        script="scripts/window_stats.py"
    output:
        "results/demo/pbs_outliers.tsv"
    shell:
        "python {input.script:q} outliers --input {input.table:q} --column pbs "
        "--top-fraction 0.01 --output {output:q}"

rule demo_rod:
    input:
        table="tests/fixtures/demo_rod.tsv",
        script="scripts/window_stats.py"
    output:
        "results/demo/rod.tsv"
    shell:
        "python {input.script:q} rod --input {input.table:q} --output {output:q}"

rule demo_fd:
    input:
        table="tests/fixtures/demo_fd.tsv",
        script="scripts/window_stats.py"
    output:
        "results/demo/fd.tsv"
    shell:
        "python {input.script:q} fd --input {input.table:q} --tail upper --sd ddof1 --output {output:q}"

rule demo_ils:
    input:
        script="scripts/window_stats.py"
    output:
        "results/demo/ils.tsv"
    shell:
        "python {input.script:q} ils --tract-bp 85400 --divergence-years 11000 "
        "--generation-years 3 --recombination-rate 1.5e-8 --branch-factor 2 --output {output:q}"

rule demo_provenance:
    input:
        outputs=expand("results/demo/{name}.tsv", name=["pbs", "pbs_outliers", "rod", "fd", "ils"]),
        script="scripts/record_provenance.py",
        config_file="config/config.test.yaml"
    output:
        "results/demo/provenance.json"
    params:
        paths=lambda wildcards, input: " ".join("--input " + shlex.quote(str(p)) for p in [*input.outputs, input.config_file])
    shell:
        "python {input.script:q} {params.paths} --output {output:q}"
