#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import unittest
import os
from modules.general import memory_usage_resource, s_encoding, sam_to_ids

OUT_MT = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data",
    "general",
    "SAMD00077852_1_MT_outmt.sam.gz"
)

class TestMemoryUsageResource(unittest.TestCase):
    def test_memory_usage_resource(self):
        # Given/When
        usage = memory_usage_resource()

        # Then
        self.assertIsInstance(usage, float)


class TestSEncoding(unittest.TestCase):
    def test_s_encoding_string(self):
        # Given
        s = "string"
        expect = "string"

        # When
        result = s_encoding(s)

        # Then
        self.assertEqual(expect, result)
        self.assertIsInstance(result, str)

    def test_s_encoding_bytes(self):
        # Given
        s = bytes("bytes", encoding="utf-8")
        expect = s.decode("utf-8")

        # When
        result = s_encoding(s)

        # Then
        self.assertEqual(expect, result)
        self.assertIsInstance(result, str)

    def test_s_encoding_invalid(self):
        # Given/When
        s = 420

        # Then
        self.assertRaises(TypeError, s_encoding, s)

class TestProcessAlignments(unittest.TestCase):
    def test_sam_to_ids(self):
        # Given/When
        r = sam_to_ids(samfile=OUT_MT, return_files=False, return_dict=True)
        # Then
        self.assertIsInstance(r, dict)
        self.assertEqual(len(r), 7060)
        self.assertEqual(len([i for i in r if r[i].bitwise_status == True]), 7060)
        self.assertEqual(len([i for i in r if r[i].paired_status == "PE"]), 4822)
        self.assertEqual(len([i for i in r if r[i].paired_status == "SE"]), 2238)

