
import unittest

from src.decay_physics import (
    remaining_nuclei,
    half_life,
    activity,
)


class TestDecayPhysics(unittest.TestCase):

    def test_initial_population(self):
        self.assertAlmostEqual(
            remaining_nuclei(1000, 0.1, 0),
            1000,
        )

    def test_half_life(self):
        t_half = half_life(0.1)
        self.assertAlmostEqual(t_half, 6.9314718, places=6)

        population = remaining_nuclei(1000, 0.1, t_half)
        self.assertAlmostEqual(population, 500, places=8)

    def test_initial_activity(self):
        self.assertAlmostEqual(activity(1000, 0.1, 0), 100)

    def test_invalid_decay_constant(self):
        with self.assertRaises(ValueError):
            half_life(0)

    def test_negative_time(self):
        with self.assertRaises(ValueError):
            remaining_nuclei(1000, 0.1, -1)


if __name__ == "__main__":
    unittest.main()
