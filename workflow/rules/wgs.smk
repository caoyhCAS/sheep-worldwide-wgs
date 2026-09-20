rule raw_fastqc:
    input:
        lambda wc: SAMPLES[wc.sample]["r1" if wc.mate == "R1" else "r2"]
    output:
        html="results/wgs/fastqc/{sample}_{mate}_fastqc.html",
        zip="results/wgs/fastqc/{sample}_{mate}_fastqc.zip"
    wildcard_constraints:
        mate="R[12]"
    threads: 1
    log: "logs/wgs/fastqc/{sample}_{mate}.log"
    conda: "../envs/reads.yaml"
    params:
        stem=lambda wc: f"{wc.sample}_{wc.mate}"
    shell:
        "python workflow/scripts/run_fastqc.py --input {input:q} --stem {params.stem:q} "
        "--html {output.html:q} --zip {output.zip:q} > {log:q} 2>&1"


rule trim_reads:
    input:
        r1=lambda wc: SAMPLES[wc.sample]["r1"],
        r2=lambda wc: SAMPLES[wc.sample]["r2"]
    output:
        r1="results/wgs/trim/{sample}.R1.paired.fq.gz",
        r2="results/wgs/trim/{sample}.R2.paired.fq.gz",
        u1=temp("results/wgs/trim/{sample}.R1.unpaired.fq.gz"),
        u2=temp("results/wgs/trim/{sample}.R2.unpaired.fq.gz")
    params: steps=TRIM_STEPS
    threads: THREADS
    log: "logs/wgs/trim/{sample}.log"
    conda: "../envs/reads.yaml"
    shell:
        "trimmomatic PE -threads {threads} {input.r1:q} {input.r2:q} "
        "{output.r1:q} {output.u1:q} {output.r2:q} {output.u2:q} {params.steps:q} > {log:q} 2>&1"


rule align_reads:
    input:
        r1="results/wgs/trim/{sample}.R1.paired.fq.gz",
        r2="results/wgs/trim/{sample}.R2.paired.fq.gz",
        ref=REFERENCE,
        bwa=[REFERENCE + suffix for suffix in [".amb", ".ann", ".bwt", ".pac", ".sa"]]
    output: temp("results/wgs/alignment/{sample}.sorted.bam")
    params:
        rg=lambda wc: r"@RG\tID:" + wc.sample + r"\tSM:" + wc.sample + r"\tPL:ILLUMINA"
    threads: THREADS
    log: "logs/wgs/align/{sample}.log"
    conda: "../envs/alignment.yaml"
    shell:
        "(bwa mem -t {threads} -R {params.rg:q} {input.ref:q} {input.r1:q} {input.r2:q} | "
        "samtools sort -@ {threads} -o {output:q} -) 2> {log:q}"


rule remove_duplicates:
    input: "results/wgs/alignment/{sample}.sorted.bam"
    output:
        bam="results/wgs/alignment/{sample}.dedup.bam",
        bai="results/wgs/alignment/{sample}.dedup.bai",
        metrics="results/wgs/alignment/{sample}.duplicate_metrics.txt"
    log: "logs/wgs/markduplicates/{sample}.log"
    conda: "../envs/gatk.yaml"
    resources: mem_mb=16000
    shell:
        "gatk --java-options '-Xmx12g' MarkDuplicates -I {input:q} -O {output.bam:q} "
        "-M {output.metrics:q} --REMOVE_DUPLICATES true --CREATE_INDEX true > {log:q} 2>&1"


rule measure_coverage:
    input:
        bam="results/wgs/alignment/{sample}.dedup.bam",
        bai="results/wgs/alignment/{sample}.dedup.bai",
        fai=REFERENCE + ".fai"
    output: "results/wgs/coverage/{sample}.json"
    log: "logs/wgs/coverage/{sample}.log"
    conda: "../envs/alignment.yaml"
    params: sample=lambda wc: wc.sample
    shell:
        "samtools depth -aa {input.bam:q} 2> {log:q} | "
        "python workflow/scripts/check_coverage.py measure --sample {params.sample:q} --fai {input.fai:q} "
        "--output {output:q} 2>> {log:q}"


rule coverage_gate:
    input: "results/wgs/coverage/{sample}.json"
    output: "results/wgs/coverage/{sample}.gt10.ok"
    conda: "../envs/python.yaml"
    shell:
        "python workflow/scripts/check_coverage.py gate --input {input:q} --threshold 10 --output {output:q}"


rule haplotype_caller:
    input:
        bam="results/wgs/alignment/{sample}.dedup.bam",
        bai="results/wgs/alignment/{sample}.dedup.bai",
        gate="results/wgs/coverage/{sample}.gt10.ok",
        ref=REFERENCE,
        fai=REFERENCE + ".fai",
        dictionary=REFERENCE_DICT
    output:
        vcf="results/wgs/gvcf/{sample}.g.vcf.gz",
        index="results/wgs/gvcf/{sample}.g.vcf.gz.tbi"
    log: "logs/wgs/haplotypecaller/{sample}.log"
    conda: "../envs/gatk.yaml"
    resources: mem_mb=16000
    shell:
        "gatk --java-options '-Xmx12g' HaplotypeCaller -R {input.ref:q} -I {input.bam:q} "
        "-ERC GVCF -O {output.vcf:q} > {log:q} 2>&1"


rule combine_gvcfs:
    input:
        vcfs=expand("results/wgs/gvcf/{sample}.g.vcf.gz", sample=CALL_SAMPLES),
        indexes=expand("results/wgs/gvcf/{sample}.g.vcf.gz.tbi", sample=CALL_SAMPLES),
        ref=REFERENCE,
        fai=REFERENCE + ".fai",
        dictionary=REFERENCE_DICT
    output:
        vcf="results/wgs/cohort.g.vcf.gz",
        index="results/wgs/cohort.g.vcf.gz.tbi"
    log: "logs/wgs/combine_gvcfs.log"
    conda: "../envs/gatk.yaml"
    resources: mem_mb=32000
    shell:
        "python workflow/scripts/combine_gvcfs.py --reference {input.ref:q} --output {output.vcf:q} "
        "{input.vcfs:q} > {log:q} 2>&1"


rule genotype_gvcfs:
    input:
        vcf="results/wgs/cohort.g.vcf.gz",
        index="results/wgs/cohort.g.vcf.gz.tbi",
        ref=REFERENCE,
        fai=REFERENCE + ".fai",
        dictionary=REFERENCE_DICT
    output:
        vcf="results/wgs/cohort.raw.vcf.gz",
        index="results/wgs/cohort.raw.vcf.gz.tbi"
    log: "logs/wgs/genotype_gvcfs.log"
    conda: "../envs/gatk.yaml"
    resources: mem_mb=32000
    shell:
        "gatk --java-options '-Xmx24g' GenotypeGVCFs -R {input.ref:q} -V {input.vcf:q} "
        "-O {output.vcf:q} > {log:q} 2>&1"


rule select_variant_type:
    input:
        vcf="results/wgs/cohort.raw.vcf.gz",
        index="results/wgs/cohort.raw.vcf.gz.tbi",
        ref=REFERENCE,
        fai=REFERENCE + ".fai",
        dictionary=REFERENCE_DICT
    output:
        vcf="results/wgs/cohort.raw.{kind}.vcf.gz",
        index="results/wgs/cohort.raw.{kind}.vcf.gz.tbi"
    wildcard_constraints: kind="snps|indels"
    params: variant_type=lambda wc: "SNP" if wc.kind == "snps" else "INDEL"
    log: "logs/wgs/select/{kind}.log"
    conda: "../envs/gatk.yaml"
    shell:
        "gatk SelectVariants -R {input.ref:q} -V {input.vcf:q} "
        "--select-type-to-include {params.variant_type:q} -O {output.vcf:q} > {log:q} 2>&1"


rule hard_filter:
    input:
        vcf="results/wgs/cohort.raw.{kind}.vcf.gz",
        index="results/wgs/cohort.raw.{kind}.vcf.gz.tbi",
        ref=REFERENCE,
        fai=REFERENCE + ".fai",
        dictionary=REFERENCE_DICT
    output:
        vcf="results/wgs/cohort.filtered.{kind}.vcf.gz",
        index="results/wgs/cohort.filtered.{kind}.vcf.gz.tbi"
    wildcard_constraints: kind="snps|indels"
    params:
        filters=lambda wc: [arg for name, expression in (
            [("QUAL30", "QUAL < 30.0"), ("QD2", "QD < 2.0"), ("MQ40", "MQ < 40.0"),
             ("FS60", "FS > 60.0"), ("SOR3", "SOR > 3.0"), ("MQRankSum", "MQRankSum < -12.5"),
             ("ReadPosRankSum", "ReadPosRankSum < -8.0")]
            if wc.kind == "snps" else [("QD2", "QD < 2.0"), ("QUAL30", "QUAL < 30.0"),
                                       ("FS200", "FS > 200.0"), ("ReadPosRankSum", "ReadPosRankSum < -20.0")])
            for arg in ["--filter-name", name, "--filter-expression", expression]],
        missing="true" if MISSING_POLICY == "fail" else "false"
    log: "logs/wgs/filter/{kind}.log"
    conda: "../envs/gatk.yaml"
    shell:
        "gatk VariantFiltration -R {input.ref:q} -V {input.vcf:q} {params.filters:q} "
        "--missing-values-evaluate-as-failing {params.missing} "
        "-O {output.vcf:q} > {log:q} 2>&1"


rule pass_variants:
    input:
        vcf="results/wgs/cohort.filtered.{kind}.vcf.gz",
        index="results/wgs/cohort.filtered.{kind}.vcf.gz.tbi",
        ref=REFERENCE,
        fai=REFERENCE + ".fai",
        dictionary=REFERENCE_DICT
    output:
        vcf="results/wgs/cohort.pass.{kind}.vcf.gz",
        index="results/wgs/cohort.pass.{kind}.vcf.gz.tbi"
    wildcard_constraints: kind="snps|indels"
    params: allelicity=lambda wc: ["--restrict-alleles-to", "BIALLELIC"] if wc.kind == "snps" else []
    log: "logs/wgs/pass/{kind}.log"
    conda: "../envs/gatk.yaml"
    shell:
        "gatk SelectVariants -R {input.ref:q} -V {input.vcf:q} --exclude-filtered true "
        "{params.allelicity:q} -O {output.vcf:q} > {log:q} 2>&1"
