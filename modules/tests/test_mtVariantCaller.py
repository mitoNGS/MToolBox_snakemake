#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste and Domenico Simone
import os
from types import SimpleNamespace
import unittest

from ..mtVariantCaller import (
    parse_mismatches_from_cigar_md,
    allele_strand_counter,
    allele_strand_updater,
    mismatch_detection,
    parse_coverage_data_file
)

r = "A00181:108:HLFMYDSXX:2:2205:3495:11303\t161\tk127_149\t569\t40\t146M1D5M\t=\t564\t-157\tCGTAACGGTTTGCTCCGTCTGACACGGCGGTTCCTTATCGAGTTGGTGTTCCCGGGCATCGTGGGCGCCGGGGGAGTTGTGAATGGCGGTACAATACCCGACGAGGAAAAATACCATGATGTTTACGCGCCATTTCATGTTGATGAGATCG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFF\tAS:i:-23\tXN:i:0\tXM:i:3\tXO:i:1\tXG:i:1\tNM:i:4\tMD:Z:32T113^C1G2T0\tYS:i:-28\tYT:Z:DP"

mismatch_dict = {601 : SimpleNamespace(DP=24, POS=601, REF='T', allele_DP=[24], allele_strand_count=[[24, 0]], alleles=['C'])}
test_sam = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data",
    "test.sam"
)
test_bam_cov = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data",
    "test.bam.cov"
) 

class TestSNPcalling(unittest.TestCase):

    def setUp(self) -> None:
        self.r = r
        self.mismatch_dict = mismatch_dict

    def test_parse_mismatches_from_cigar_md_default(self):
        # Given
        expected = ([601], [32], ['T'], ['C'], [37], '+')
        
        # When
        result = parse_mismatches_from_cigar_md(self.r)
        
        # Then
        self.assertEqual(expected, result)

    def test_parse_mismatches_from_cigar_md_tail_zero(self):
        # Given
        expected = ([601, 717, 720], [32, 148, 151], ['T', 'G', 'T'], ['C', 'A', 'G'], [37, 37, 37], '+')
        
        # When
        result = parse_mismatches_from_cigar_md(self.r, tail_mismatch=0)
        
        # Then
        self.assertEqual(expected, result)
    
    def test_allele_strand_counter(self):
        # Given
        expected = [[1, 0], [0, 1]]
        
        # When
        result = [allele_strand_counter("+"), allele_strand_counter("-")]
        
        # Then
        self.assertEqual(expected, result)

    def test_allele_strand_updater(self):
        # Given
        expected = [25, 5]
        # When
        result = allele_strand_updater([0, 1], [25, 4])
        # Then
        self.assertEqual(expected, result)

    def test_mismatch_detection(self):
        # Given
        expected = self.mismatch_dict
        # When
        result = mismatch_detection(sam=test_sam, coverage_data=parse_coverage_data_file(test_bam_cov))
        # Then
        self.assertEqual(expected, result)
        
    # 
    # def test_get_mt_genomes(self):
    #     # Given
    #     expected = ["NC_001323.1"]
    # 
    #     # When
    #     result = get_mt_genomes(self.analysis_tab)
    # 
    #     # Then
    #     self.assertEqual(expected, result)

    # def test_fastqc_outputs_raw(self):
    #     # Given
    #     expected = [el.format("results/fastqc_raw")
    #                 for el in FASTQC_OUT]
    # 
    #     # When
    #     result = fastqc_outputs(datasets_tab=self.datasets_tab,
    #                             analysis_tab=self.analysis_tab,
    #                             out="raw")
    # 
    #     # Then
    #     self.assertEqual(expected, result)
    # 
    # def test_fastqc_outputs_filtered(self):
    #     # Given
    #     expected = [el.format("results/fastqc_filtered")
    #                 for el in FASTQC_OUT_FILT]
    # 
    #     # When
    #     result = fastqc_outputs(datasets_tab=self.datasets_tab,
    #                             analysis_tab=self.analysis_tab,
    #                             out="filtered")
    # 
    #     # Then
    #     self.assertEqual(sorted(expected), sorted(result))
    # 
    # def test_fastqc_outputs_error(self):
    #     # Given/When
    #     with self.assertRaises(ValueError):
    #         fastqc_outputs(datasets_tab=self.datasets_tab,
    #                        analysis_tab=self.analysis_tab,
    #                        out="test")
    # 
    # def test_get_genome_vcf_files(self):
    #     # Given
    #     expected = ["results/vcf/NC_001323.1_GCF_000002315.5.vcf"]
    # 
    #     # When
    #     result = get_genome_vcf_files(df=self.analysis_tab)
    # 
    #     # Then
    #     self.assertEqual(expected, result)
    # 
    # def test_get_bed_files(self):
    #     # Given
    #     expected = [
    #         "results/5517_hypo/5517_hypo_NC_001323.1_GCF_000002315.5.bed",
    #         "results/5517_liver/5517_liver_NC_001323.1_GCF_000002315.5.bed"
    #     ]
    # 
    #     # When
    #     result = get_bed_files(df=self.analysis_tab)
    # 
    #     # Then
    #     self.assertEqual(expected, result)
    # 
    # def test_get_fasta_files(self):
    #     # Given
    #     expected = [
    #         "results/5517_hypo/5517_hypo_NC_001323.1_GCF_000002315.5.fasta",
    #         "results/5517_liver/5517_liver_NC_001323.1_GCF_000002315.5.fasta"
    #     ]
    # 
    #     # When
    #     result = get_fasta_files(df=self.analysis_tab)
    # 
    #     # Then
    #     self.assertEqual(expected, result)
    # 
    # def test_get_genome_files(self):
    #     # Given
    #     expected = ["GCF_000002315.5_GRCg6a_genomic_mt.fna"]
    # 
    #     # When
    #     result = get_genome_files(df=self.reference_tab,
    #                               ref_genome_mt="NC_001323.1",
    #                               field="ref_genome_mt_file")
    # 
    #     # Then
    #     self.assertEqual(expected, result)
