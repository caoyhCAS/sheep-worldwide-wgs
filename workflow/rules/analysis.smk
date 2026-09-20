rule prepare_autosomal_vcf:
    input:
        vcf=SOURCE_VCF,
        index=SOURCE_VCF + ".tbi",
        rename=CONTIG_MAP
    output:
        vcf="results/analysis/autosomal.pass.snps.vcf.gz",
        index="results/analysis/autosomal.pass.snps.vcf.gz.tbi"
    params: autosomes=AUTOSOMES
    log: "logs/analysis/prepare_vcf.log"
    conda: "../envs/variants.yaml"
    shell:
        "(bcftools annotate --rename-chrs {input.rename:q} --set-id '%CHROM:%POS:%REF:%FIRST_ALT' "
        "-Ou {input.vcf:q} | bcftools view --apply-filters PASS --types snps --min-alleles 2 --max-alleles 2 "
        "--targets {params.autosomes:q} -Oz -o {output.vcf:q}) 2> {log:q}\n"
        "bcftools index --tbi {output.vcf:q} 2>> {log:q}"


rule validate_analysis_cohort:
    input:
        vcf="results/analysis/autosomal.pass.snps.vcf.gz",
        manifest=ANALYSIS["populations"],
        keeps=[row["keep"] for row in POPULATIONS.values()],
        unrelated=[UNRELATED_KEEP] if "population" in MODULES else [],
        ld=[row["ld_keep"] for row in POPULATIONS.values()] if "population" in MODULES and LD_ENABLED else []
    output:
        report="results/analysis/cohort_validation.json"
    params:
        unrelated=["--unrelated", UNRELATED_KEEP] if "population" in MODULES else [],
        ld=["--check-ld"] if "population" in MODULES and LD_ENABLED else []
    conda: "../envs/python.yaml"
    shell:
        "python workflow/scripts/validate_cohort.py --vcf {input.vcf:q} --populations {input.manifest:q} "
        "{params.unrelated:q} {params.ld:q} --output {output.report:q}"
