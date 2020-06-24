import snakemake, os, sys
import mtoolnote
from modules.config_parsers import (
    get_analysis_species, parse_config_tabs
)

analysis_tab, reference_tab, datasets_tab = parse_config_tabs(analysis_tab_file="data/analysis.tab", reference_tab_file="data/reference_genomes.tab", datasets_tab_file="data/datasets.tab")
species = get_analysis_species(ref_genome_mt, reference_tab=None, config_species=None)

### config files and configuration
configfile: "config.yaml"

include: "variant_calling.snakefile"

rule all_annotation:
    input: get_genome_vcf_files(analysis_tab, annotation=True)

rule annotate_vcf:
    input:
        vcf = "results/vcf/{ref_genome_mt}_{ref_genome_n}.vcf"
    output:
        vcf = "results/vcf/{ref_genome_mt}_{ref_genome_n}.annotated.vcf"
    params:
        species = species
    run:
        mtoolnote.annotate(input.vcf, output.vcf, species=params.species)