rule plink_import:
    input:
        vcf="results/analysis/autosomal.pass.snps.vcf.gz",
        validation="results/analysis/cohort_validation.json"
    output:
        bed="results/population/all.bed",
        bim="results/population/all.bim",
        fam="results/population/all.fam"
    params: prefix="results/population/all"
    log: "logs/population/plink_import.log"
    conda: "../envs/plink.yaml"
    shell:
        "plink --vcf {input.vcf:q} --double-id --chr-set 26 no-xy --keep-allele-order "
        "--make-bed --out {params.prefix:q} > {log:q} 2>&1"


rule plink_unrelated:
    input:
        bed="results/population/all.bed",
        bim="results/population/all.bim",
        fam="results/population/all.fam",
        keep=UNRELATED_KEEP
    output:
        bed="results/population/unrelated.bed",
        bim="results/population/unrelated.bim",
        fam="results/population/unrelated.fam"
    params: source="results/population/all", prefix="results/population/unrelated"
    log: "logs/population/unrelated.log"
    conda: "../envs/plink.yaml"
    shell:
        "plink --bfile {params.source:q} --chr-set 26 no-xy --keep {input.keep:q} "
        "--keep-allele-order --make-bed --out {params.prefix:q} > {log:q} 2>&1"


rule population_qc:
    input:
        bed="results/population/unrelated.bed",
        bim="results/population/unrelated.bim",
        fam="results/population/unrelated.fam"
    output:
        bed="results/population/qc.bed",
        bim="results/population/qc.bim",
        fam="results/population/qc.fam"
    params: source="results/population/unrelated", prefix="results/population/qc"
    log: "logs/population/qc.log"
    conda: "../envs/plink.yaml"
    shell:
        "plink --bfile {params.source:q} --chr-set 26 no-xy --maf 0.05 --geno 0.1 --hwe 0.001 "
        "--keep-allele-order --make-bed --out {params.prefix:q} > {log:q} 2>&1"


rule ld_prune_list:
    input:
        bed="results/population/qc.bed",
        bim="results/population/qc.bim",
        fam="results/population/qc.fam"
    output:
        keep="results/population/ld.prune.in",
        removed="results/population/ld.prune.out"
    params: source="results/population/qc", prefix="results/population/ld"
    log: "logs/population/prune.log"
    conda: "../envs/plink.yaml"
    shell:
        "plink --bfile {params.source:q} --chr-set 26 no-xy --indep-pairwise 50 5 0.2 "
        "--out {params.prefix:q} > {log:q} 2>&1"


rule apply_ld_pruning:
    input:
        bed="results/population/qc.bed",
        bim="results/population/qc.bim",
        fam="results/population/qc.fam",
        keep="results/population/ld.prune.in"
    output:
        bed="results/population/pruned.bed",
        bim="results/population/pruned.bim",
        fam="results/population/pruned.fam"
    params: source="results/population/qc", prefix="results/population/pruned"
    log: "logs/population/apply_pruning.log"
    conda: "../envs/plink.yaml"
    shell:
        "plink --bfile {params.source:q} --chr-set 26 no-xy --extract {input.keep:q} "
        "--keep-allele-order --make-bed --out {params.prefix:q} > {log:q} 2>&1"


rule ibs_distance:
    input:
        bed="results/population/pruned.bed",
        bim="results/population/pruned.bim",
        fam="results/population/pruned.fam"
    output:
        matrix="results/population/ibs.mdist",
        ids="results/population/ibs.mdist.id"
    params: source="results/population/pruned", prefix="results/population/ibs"
    log: "logs/population/ibs.log"
    conda: "../envs/plink.yaml"
    shell:
        "plink --bfile {params.source:q} --chr-set 26 no-xy --distance square 1-ibs "
        "--out {params.prefix:q} > {log:q} 2>&1"


rule population_heterozygosity:
    input:
        bed="results/population/unrelated.bed",
        bim="results/population/unrelated.bim",
        fam="results/population/unrelated.fam"
    output: "results/population/heterozygosity.het"
    params: source="results/population/unrelated", prefix="results/population/heterozygosity"
    log: "logs/population/heterozygosity.log"
    conda: "../envs/plink.yaml"
    shell:
        "plink --bfile {params.source:q} --chr-set 26 no-xy --het "
        "--out {params.prefix:q} > {log:q} 2>&1"


rule individual_pi:
    input:
        vcf="results/analysis/autosomal.pass.snps.vcf.gz",
        validation="results/analysis/cohort_validation.json"
    output: "results/population/individual_pi/{individual}.windowed.pi"
    params: prefix=lambda wc: f"results/population/individual_pi/{wc.individual}"
    log: "logs/population/individual_pi/{individual}.log"
    conda: "../envs/variants.yaml"
    shell:
        "vcftools --gzvcf {input.vcf:q} --indv {wildcards.individual:q} --window-pi 1000000 "
        "--window-pi-step 1000000 --out {params.prefix:q} > {log:q} 2>&1"


rule population_pi:
    input:
        vcf="results/analysis/autosomal.pass.snps.vcf.gz",
        validation="results/analysis/cohort_validation.json",
        keep=lambda wc: POPULATIONS[wc.population]["keep"]
    output: "results/population/pi/{population}.windowed.pi"
    params: prefix=lambda wc: f"results/population/pi/{wc.population}"
    log: "logs/population/pi/{population}.log"
    conda: "../envs/variants.yaml"
    shell:
        "vcftools --gzvcf {input.vcf:q} --keep {input.keep:q} --window-pi 1000000 "
        "--window-pi-step 1000000 --out {params.prefix:q} > {log:q} 2>&1"


if ROH_ENABLED:
    rule population_roh:
        input:
            bed="results/population/unrelated.bed",
            bim="results/population/unrelated.bim",
            fam="results/population/unrelated.fam"
        output:
            segments="results/population/roh.hom",
            individuals="results/population/roh.hom.indiv"
        params:
            source="results/population/unrelated", prefix="results/population/roh", extra=ROH_ARGS
        log: "logs/population/roh.log"
        conda: "../envs/plink.yaml"
        shell:
            "plink --bfile {params.source:q} --chr-set 26 no-xy --homozyg --homozyg-density 50 "
            "--homozyg-window-het 1 --homozyg-window-snp 50 {params.extra:q} "
            "--out {params.prefix:q} > {log:q} 2>&1"


if LD_ENABLED:
    rule population_ld:
        input:
            vcf="results/analysis/autosomal.pass.snps.vcf.gz",
            validation="results/analysis/cohort_validation.json",
            keep=lambda wc: POPULATIONS[wc.population]["ld_keep"]
        output: "results/population/ld/{population}.stat.gz"
        params: prefix=lambda wc: f"results/population/ld/{wc.population}"
        log: "logs/population/ld/{population}.log"
        conda: "../envs/poplddecay.yaml"
        shell:
            "PopLDdecay -InVCF {input.vcf:q} -SubPop {input.keep:q} -MaxDist 300 "
            "-OutStat {params.prefix:q} > {log:q} 2>&1"
