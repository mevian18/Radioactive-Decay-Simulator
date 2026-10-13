
"""Monte Carlo simulation of radioactive decay."""

import numpy as np


def simulate_decay(
    initial_nuclei,
    decay_constant,
    time_step,
    total_time,
    seed=None,
):
    """Simulate radioactive decay using random sampling."""

    # Validate the simulation inputs
    if (
        not isinstance(initial_nuclei, (int, np.integer))
        or isinstance(initial_nuclei, (bool, np.bool_))
    ):
        raise ValueError("Initial nuclei must be a non-negative integer.")

    if initial_nuclei < 0:
        raise ValueError("Initial nuclei cannot be negative.")

    if not np.isfinite([decay_constant, time_step, total_time]).all():
        raise ValueError("Simulation parameters must be finite.")

    if decay_constant < 0:
        raise ValueError("Decay constant cannot be negative.")

    if time_step <= 0:
        raise ValueError("Time step must be positive.")

    if total_time < 0:
        raise ValueError("Total time cannot be negative.")

    # Create the random number generator
    rng = np.random.default_rng(seed)

    # Generate time points, including the exact final time
    times = np.arange(0, total_time, time_step)
    times = np.append(times, total_time)

    # Store the number of surviving nuclei at each time
    populations = np.zeros(len(times), dtype=int)
    populations[0] = initial_nuclei

    # Simulate decay over each time interval
    for step in range(len(times) - 1):
        current_population = populations[step]

        actual_time_step = times[step + 1] - times[step]

        probability_of_decay = 1 - np.exp(
            -decay_constant * actual_time_step
        )

        decayed = rng.binomial(
            current_population,
            probability_of_decay,
        )

        populations[step + 1] = current_population - decayed

    return times, populations
