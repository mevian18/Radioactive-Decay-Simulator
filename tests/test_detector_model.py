
"""Tests for detector counting and activity estimation."""
from statistics import NormalDist

import numpy as np
import unittest

from src.detector_model import (
    estimate_activity,
    theoretical_activity_std,
    critical_level_counts,
    detection_limit_counts,
    minimum_detectable_activity,
)


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
def test_currie_identity():
    # for alpha == beta: L_D = k^2 + 2 * L_C
    k = NormalDist().inv_cdf(0.95)
    lc = critical_level_counts(5.0, 60.0, 12.0)
    ld = detection_limit_counts(5.0, 60.0, 12.0)
    assert abs(ld - (k**2 + 2 * lc)) < 1e-9


def test_detection_limit_monte_carlo():
    rng = np.random.default_rng(2024)
    eps, rb, ts, tb, n = 0.8, 5.0, 60.0, 60.0, 100_000
    lc = critical_level_counts(rb, ts, tb)
    mda = minimum_detectable_activity(eps, rb, ts, tb)
    g0 = rng.poisson(rb * ts, n)
    b0 = rng.poisson(rb * tb, n)
    g1 = rng.poisson((eps * mda + rb) * ts, n)
    b1 = rng.poisson(rb * tb, n)
    false_positive = np.mean(g0 - (ts / tb) * b0 > lc)
    detection = np.mean(g1 - (ts / tb) * b1 > lc)
    assert abs(false_positive - 0.05) < 0.006
    assert abs(detection - 0.95) < 0.006