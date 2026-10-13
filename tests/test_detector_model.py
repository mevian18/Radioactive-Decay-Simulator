
"""Tests for detector counting and activity estimation."""

import unittest

from src.detector_model import estimate_activity


class TestDetectorModel(unittest.TestCase):

    def test_known_activity(self):
        result = estimate_activity(
            gross_counts=850,
            background_counts=50,
            source_time=10,
            background_time=10,
            efficiency=0.8,
        )

        # Net rate = 85 - 5 = 80 counts/s
        # Activity estimate = 80 / 0.8 = 100 Bq
        self.assertAlmostEqual(
            result["estimated_activity_bq"], 100.0
        )

    def test_zero_net_counts(self):
        result = estimate_activity(
            gross_counts=50,
            background_counts=50,
            source_time=10,
            background_time=10,
            efficiency=0.5,
        )

        self.assertAlmostEqual(
            result["estimated_activity_bq"], 0.0
        )

    def test_negative_net_estimate_is_preserved(self):
        result = estimate_activity(
            gross_counts=40,
            background_counts=50,
            source_time=10,
            background_time=10,
            efficiency=0.5,
        )

        self.assertLess(
            result["estimated_activity_bq"], 0.0
        )

    def test_invalid_efficiency(self):
        with self.assertRaises(ValueError):
            estimate_activity(100, 10, 10, 10, 0)

    def test_invalid_measurement_time(self):
        with self.assertRaises(ValueError):
            estimate_activity(100, 10, 0, 10, 0.5)

    def test_negative_counts(self):
        with self.assertRaises(ValueError):
            estimate_activity(-1, 10, 10, 10, 0.5)


if __name__ == "__main__":
    unittest.main()
