#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import os
import unittest

import pandas as pd

from ..config_parsers import fastqc_raw_outputs, get_mt_genomes

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
FASTQC_OUT = [
    "results/fastqc_raw/5517_hypo_1.R1_fastqc.html", 
    "results/fastqc_raw/5517_hypo_1.R2_fastqc.html", 
    "results/fastqc_raw/5517_hypo_2.R1_fastqc.html", 
    "results/fastqc_raw/5517_hypo_2.R2_fastqc.html", 
    "results/fastqc_raw/5517_liver_1.R1_fastqc.html", 
    "results/fastqc_raw/5517_liver_1.R2_fastqc.html", 
    "results/fastqc_raw/5517_liver_2.R1_fastqc.html", 
    "results/fastqc_raw/5517_liver_2.R2_fastqc.html"
]


class TestConfigParsers(unittest.TestCase):

    def setUp(self) -> None:
        self.analysis_tab = pd.read_table(ANALYSIS_TAB, sep="\t",
                                          comment="#")
        self.datasets_tab = pd.read_table(DATASETS_TAB, sep="\t",
                                          comment="#")

    def test_get_mt_genomes(self):
        # Given
        expected = ["NC_001323.1"]

        # When
        result = get_mt_genomes(self.analysis_tab)

        # Then
        self.assertEqual(expected, result)

    def test_fastqc_raw_outputs(self):
        # Given/When
        result = fastqc_raw_outputs(datasets_tab=self.datasets_tab,
                                    analysis_tab=self.analysis_tab)

        # Then
        self.assertEqual(FASTQC_OUT, result)
