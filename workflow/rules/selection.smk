if GLOBAL_FST:
    rule selection_global_fst:
        input:
            vcf="results/analysis/autosomal.pass.snps.vcf.gz",
            validation="results/analysis/cohort_validation.json",
            keeps=[row["keep"] for row in POPULATIONS.values()]
        output: "results/selection/global_fst.windowed.weir.fst"
        params:
            prefix="results/selection/global_fst",
            window=FST_WINDOW,
            step=FST_STEP,
            populations=lambda wc, input: " ".join("--weir-fst-pop " + shlex.quote(str(path)) for path in input.keeps)
        log: "logs/selection/global_fst.log"
        conda: "../envs/variants.yaml"
        shell:
            "vcftools --gzvcf {input.vcf:q} {params.populations} "
            "--fst-window-size {params.window} --fst-window-step {params.step} "
            "--out {params.prefix:q} > {log:q} 2>&1"


rule selection_fst:
    input:
        vcf="results/analysis/autosomal.pass.snps.vcf.gz",
        validation="results/analysis/cohort_validation.json",
        a=lambda wc: POPULATIONS[COMPARISONS[wc.comparison][0]]["keep"],
        b=lambda wc: POPULATIONS[COMPARISONS[wc.comparison][1]]["keep"]
    output: "results/selection/fst/{comparison}.windowed.weir.fst"
    params:
        prefix=lambda wc: f"results/selection/fst/{wc.comparison}",
        window=FST_WINDOW,
        step=FST_STEP
    log: "logs/selection/fst/{comparison}.log"
    conda: "../envs/variants.yaml"
    shell:
        "vcftools --gzvcf {input.vcf:q} --weir-fst-pop {input.a:q} --weir-fst-pop {input.b:q} "
        "--fst-window-size {params.window} --fst-window-step {params.step} "
        "--out {params.prefix:q} > {log:q} 2>&1"


rule selection_pi:
    input:
        vcf="results/analysis/autosomal.pass.snps.vcf.gz",
        validation="results/analysis/cohort_validation.json",
        keep=lambda wc: POPULATIONS[wc.population]["keep"]
    output: "results/selection/pi/{population}.windowed.pi"
    params: prefix=lambda wc: f"results/selection/pi/{wc.population}"
    log: "logs/selection/pi/{population}.log"
    conda: "../envs/variants.yaml"
    shell:
        "vcftools --gzvcf {input.vcf:q} --keep {input.keep:q} --window-pi 10000 --window-pi-step 10000 "
        "--out {params.prefix:q} > {log:q} 2>&1"
