#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import unittest

from scripts.filter_HF import to_list


class TestToList(unittest.TestCase):
    def test_to_list_string(self):
        # Given
        obj = "string"
        expect = ["string"]

        # When
        result = to_list(obj)

        # Then
        self.assertEqual(expect, result)
        self.assertIsInstance(result, list)

    def test_to_list_other(self):
        # Given
        obj = ["string"]

        # When
        result = to_list(obj)

        # Then
        self.assertEqual(obj, result)
