#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import os
import unittest

import pandas as pd

from ..config_parsers import (
    fastqc_outputs, get_bed_files, get_genome_files, get_fasta_files,
    get_genome_vcf_files, get_mt_genomes, get_analysis_species,
    ref_genome_mt_to_species

)

ANALYSIS_TAB = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data",
    "analysis.tab"
)
DATASETS_TAB = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data",
    "datasets.tab"
)
REFERENCE_TAB = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data",
    "reference_genomes.tab"
)
FASTQC_OUT = [
    "{}/5517_hypo_1.R1_fastqc.html",
    "{}/5517_hypo_1.R2_fastqc.html",
    "{}/5517_hypo_2.R1_fastqc.html",
    "{}/5517_hypo_2.R2_fastqc.html",
    "{}/5517_liver_1.R1_fastqc.html",
    "{}/5517_liver_1.R2_fastqc.html",
    "{}/5517_liver_2.R1_fastqc.html",
    "{}/5517_liver_2.R2_fastqc.html"
]

FASTQC_OUT_FILT = [
    "{}/5517_hypo_1_qc_R1_fastqc.html",
    "{}/5517_hypo_1_qc_R2_fastqc.html",
    "{}/5517_hypo_2_qc_R1_fastqc.html",
    "{}/5517_hypo_2_qc_R2_fastqc.html",
    "{}/5517_liver_1_qc_R1_fastqc.html",
    "{}/5517_liver_1_qc_R2_fastqc.html",
    "{}/5517_liver_2_qc_R1_fastqc.html",
    "{}/5517_liver_2_qc_R2_fastqc.html",
    "{}/5517_hypo_1_qc_U_fastqc.html",
    "{}/5517_hypo_2_qc_U_fastqc.html",
    "{}/5517_liver_1_qc_U_fastqc.html",
    "{}/5517_liver_2_qc_U_fastqc.html"
]
# FASTQC_OUT_FILT = FASTQC_OUT + [
#     "{}/5517_hypo_1_qc_U_fastqc.html",
#     "{}/5517_hypo_2_qc_U_fastqc.html",
#     "{}/5517_liver_1_qc_U_fastqc.html",
#     "{}/5517_liver_2_qc_U_fastqc.html",
# ]


class TestConfigParsers(unittest.TestCase):

    def setUp(self) -> None:
        self.analysis_tab = pd.read_table(ANALYSIS_TAB, sep="\t", comment="#")
        self.datasets_tab = pd.read_table(DATASETS_TAB, sep="\t", comment="#")
        self.reference_tab = (pd.read_table(REFERENCE_TAB, sep="\t", comment='#')
                              .set_index("ref_genome_mt", drop=False))

    def test_ref_genome_mt_to_species(self):
        # Given
        expected = "ggallus"
        # When
        result = ref_genome_mt_to_species(ref_genome_mt="NC_001323.1", reference_tab=self.reference_tab)
        # Then
        self.assertEqual(expected, result)

    def test_get_analysis_species(self):
        # Given
        expected = ["ggallus", "human"]
        # When
        result1 = get_analysis_species("NC_001323.1", reference_tab=self.reference_tab, config_species=None)
        result2 = get_analysis_species("NC_001323.1", reference_tab=self.reference_tab, config_species="human")
        # Then
        self.assertEqual(expected, [result1, result2])
        
    def test_get_mt_genomes(self):
        # Given
        expected = ["NC_001323.1"]

        # When
        result = get_mt_genomes(self.analysis_tab)

        # Then
        self.assertEqual(expected, result)

    def test_fastqc_outputs_raw(self):
        # Given
        expected = [el.format("results/fastqc_raw")
                    for el in FASTQC_OUT]

        # When
        result = fastqc_outputs(datasets_tab=self.datasets_tab,
                                analysis_tab=self.analysis_tab,
                                out="raw")

        # Then
        self.assertEqual(expected, result)

    def test_fastqc_outputs_filtered(self):
        # Given
        expected = [el.format("results/fastqc_filtered")
                    for el in FASTQC_OUT_FILT]

        # When
        result = fastqc_outputs(datasets_tab=self.datasets_tab,
                                analysis_tab=self.analysis_tab,
                                out="filtered")

        # Then
        self.assertEqual(sorted(expected), sorted(result))

    def test_fastqc_outputs_error(self):
        # Given/When
        with self.assertRaises(ValueError):
            fastqc_outputs(datasets_tab=self.datasets_tab,
                           analysis_tab=self.analysis_tab,
                           out="test")

    def test_get_genome_vcf_files(self):
        # Given
        expected = ["results/vcf/NC_001323.1_GCF_000002315.5.vcf"]

        # When
        result = get_genome_vcf_files(df=self.analysis_tab)

        # Then
        self.assertEqual(expected, result)

    def test_get_bed_files(self):
        # Given
        expected = [
            "results/5517_hypo/5517_hypo_NC_001323.1_GCF_000002315.5.bed",
            "results/5517_liver/5517_liver_NC_001323.1_GCF_000002315.5.bed"
        ]

        # When
        result = get_bed_files(df=self.analysis_tab)

        # Then
        self.assertEqual(expected, result)

    def test_get_fasta_files(self):
        # Given
        expected = [
            "results/5517_hypo/5517_hypo_NC_001323.1_GCF_000002315.5.fasta",
            "results/5517_liver/5517_liver_NC_001323.1_GCF_000002315.5.fasta"
        ]

        # When
        result = get_fasta_files(df=self.analysis_tab)

        # Then
        self.assertEqual(expected, result)

    def test_get_genome_files(self):
        # Given
        expected = ["GCF_000002315.5_GRCg6a_genomic_mt.fna"]

        # When
        result = get_genome_files(df=self.reference_tab,
                                  ref_genome_mt="NC_001323.1",
                                  field="ref_genome_mt_file")

        # Then
        self.assertEqual(expected, result)
