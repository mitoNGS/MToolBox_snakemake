#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import unittest

from modules.general import memory_usage_resource, s_encoding


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
