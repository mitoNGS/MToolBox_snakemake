import snakemake, os, sys
import mtoolnote

### config files and configuration
configfile: "config.yaml"
#localrules: all, filterSeqs

include: "variant_calling.snakefile"

rule all_annotation:
    input: get_genome_vcf_files(analysis_tab, annotation=True)

rule annotate_vcf:
    input:
        vcf = "results/vcf/{ref_genome_mt}_{ref_genome_n}.vcf"
    output:
        vcf = "results/vcf/{ref_genome_mt}_{ref_genome_n}.annotated.vcf"
    run:
        mtoolnote.annotate(input.vcf, output.vcf)