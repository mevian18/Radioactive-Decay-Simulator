
"""Analytical calculations for radioactive decay."""

import math


def remaining_nuclei(initial_nuclei, decay_constant, time):
    """Calculate the expected number of nuclei remaining.

    decay_constant is expressed in inverse units of time.
    time must use the corresponding time unit.
    """
    if initial_nuclei < 0:
        raise ValueError("Initial population cannot be negative.")
    if decay_constant < 0:
        raise ValueError("Decay constant cannot be negative.")
    if time < 0:
        raise ValueError("Time cannot be negative.")

    return initial_nuclei * math.exp(-decay_constant * time)


def half_life(decay_constant):
    """Calculate half-life in the time unit used by the decay constant."""
    if decay_constant <= 0:
        raise ValueError("Decay constant must be positive.")

    return math.log(2) / decay_constant


def activity(initial_nuclei, decay_constant, time):
    """Calculate expected activity in Bq when time is in seconds."""
    if decay_constant < 0:
        raise ValueError("Decay constant cannot be negative.")

    nuclei = remaining_nuclei(initial_nuclei, decay_constant, time)
    return decay_constant * nuclei
