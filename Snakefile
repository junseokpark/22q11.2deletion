import os

SAMPLES, = glob_wildcards("source/MSSM_{sample}_NeuN_pl.cram")
TOOLS = ["cnvnator", "gatk", "delly", "lumpy", "manta", "canvas", "meerkat"]

rule all:
    input:
        expand("results/cnvnator/{sample}.cnvnator", sample=SAMPLES),
        "results/gatk/gcnv_cohort_calls.vcf",
        expand("results/delly/{sample}.vcf", sample=SAMPLES),
        expand("results/lumpy/{sample}.vcf", sample=SAMPLES),
        expand("results/manta/{sample}/results/variants/diploidSV.vcf", sample=SAMPLES),
        expand("results/canvas/{sample}.cnv.vcf", sample=SAMPLES),
        expand("results/meerkat/{sample}/final_fusion.txt", sample=SAMPLES),
        "results/summary/variant_overlap_table.tsv"

rule cnvnator:
    input:
        cram="source/MSSM_{sample}_NeuN_pl.cram",
        ref=config["reference"]
    output:
        "results/cnvnator/{sample}.cnvnator"
    shell:
        "bash scripts/run_cnvnator.sh {input.cram} {output} {config[reference]} {config[cnvnator_data_dir]}"

rule gatk:
    input:
        crams=expand("source/MSSM_{sample}_NeuN_pl.cram", sample=SAMPLES),
        ref=config["reference"],
        ploidy_dir=config["gatk_ploidy_dir"]
    output:
        "results/gatk/gcnv_cohort_calls.vcf"
    shell:
        "bash scripts/run_gatk.sh {config[reference]} {config[gatk_ploidy_dir]} results/gatk"

rule delly:
    input:
        cram="source/MSSM_{sample}_NeuN_pl.cram",
        ref=config["reference"]
    output:
        "results/delly/{sample}.vcf"
    shell:
        "bash scripts/run_delly.sh {input.cram} {output} {input.ref}"

rule lumpy:
    input:
        cram="source/MSSM_{sample}_NeuN_pl.cram"
    output:
        "results/lumpy/{sample}.vcf"
    shell:
        "bash scripts/run_lumpy.sh {input.cram} {output}"

rule manta:
    input:
        cram="source/MSSM_{sample}_NeuN_pl.cram",
        ref=config["reference"]
    output:
        "results/manta/{sample}/results/variants/diploidSV.vcf"
    shell:
        "bash scripts/run_manta.sh {input.cram} {output} {input.ref} {wildcards.sample} {config[threads]}"

rule canvas:
    input:
        cram="source/MSSM_{sample}_NeuN_pl.cram",
        ref=config["reference"]
    output:
        "results/canvas/{sample}.cnv.vcf"
    shell:
        "bash scripts/run_canvas.sh {input.cram} {output} {input.ref} {wildcards.sample}"

rule meerkat:
    input:
        cram="source/MSSM_{sample}_NeuN_pl.cram",
        ref=config["reference"]
    output:
        "results/meerkat/{sample}/final_fusion.txt"
    shell:
        "bash scripts/run_meerkat.sh {input.cram} {output} {input.ref} {wildcards.sample}"

rule compare_variants:
    input:
        cnvnator=expand("results/cnvnator/{sample}.cnvnator", sample=SAMPLES),
        gatk="results/gatk/gcnv_cohort_calls.vcf",
        delly=expand("results/delly/{sample}.vcf", sample=SAMPLES),
        lumpy=expand("results/lumpy/{sample}.vcf", sample=SAMPLES),
        manta=expand("results/manta/{sample}/results/variants/diploidSV.vcf", sample=SAMPLES),
        canvas=expand("results/canvas/{sample}.cnv.vcf", sample=SAMPLES),
        meerkat=expand("results/meerkat/{sample}/final_fusion.txt", sample=SAMPLES)
    output:
        "results/summary/variant_overlap_table.tsv"
    shell:
        "python scripts/compare_variants.py "
        "--cnvnator {" ".join(input.cnvnator)} "
        "--gatk {input.gatk} "
        "--delly {" ".join(input.delly)} "
        "--lumpy {" ".join(input.lumpy)} "
        "--manta {" ".join(input.manta)} "
        "--canvas {" ".join(input.canvas)} "
        "--meerkat {" ".join(input.meerkat)} "
        "--output {output}"
