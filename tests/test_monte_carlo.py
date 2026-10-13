
import unittest
import numpy as np

from src.monte_carlo import simulate_decay


class TestMonteCarlo(unittest.TestCase):

    def test_initial_population(self):
        times, populations = simulate_decay(
            1000, 0.1, 0.1, 10, seed=42
        )
        self.assertEqual(populations[0], 1000)

    def test_population_never_increases(self):
        _, populations = simulate_decay(
            1000, 0.1, 0.1, 10, seed=42
        )
        self.assertTrue(np.all(np.diff(populations) <= 0))

    def test_population_is_nonnegative(self):
        _, populations = simulate_decay(
            1000, 0.1, 0.1, 10, seed=42
        )
        self.assertTrue(np.all(populations >= 0))

    def test_zero_decay_constant(self):
        _, populations = simulate_decay(
            1000, 0.0, 0.1, 10, seed=42
        )
        self.assertTrue(np.all(populations == 1000))

    def test_invalid_time_step(self):
        with self.assertRaises(ValueError):
            simulate_decay(1000, 0.1, 0, 10)

    def test_reproducible_seed(self):
        _, first = simulate_decay(1000, 0.1, 0.1, 10, seed=42)
        _, second = simulate_decay(1000, 0.1, 0.1, 10, seed=42)
        np.testing.assert_array_equal(first, second)
    
    def test_partial_final_time_step(self):
        times, populations = simulate_decay(
            1000, 0.1, 0.3, 1.0, seed=42
        )

        self.assertAlmostEqual(times[-1], 1.0)
        self.assertEqual(len(times), len(populations))
        self.assertTrue(np.all(np.diff(populations) <= 0))


if __name__ == "__main__":
    unittest.main()
