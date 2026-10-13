
"""Compare analytical and Monte Carlo radioactive decay."""

import numpy as np
import matplotlib.pyplot as plt

from src.decay_physics import remaining_nuclei
from src.monte_carlo import simulate_decay


def main():
    initial_nuclei = 1000
    decay_constant = 0.1
    time_step = 0.1
    total_time = 30.0

    # Analytical prediction
    analytical_times = np.linspace(0, total_time, 301)
    analytical_population = initial_nuclei * np.exp(
        -decay_constant * analytical_times
    )

    # Monte Carlo simulation
    mc_times, mc_population = simulate_decay(
        initial_nuclei,
        decay_constant,
        time_step,
        total_time,
        seed=42,
    )

    # Plot both results
    plt.figure(figsize=(10, 6))

    plt.plot(
        analytical_times,
        analytical_population,
        label="Analytical prediction",
        linewidth=2.5,
    )

    plt.step(
        mc_times,
        mc_population,
        where="post",
        label="Monte Carlo simulation",
    )

    plt.xlabel("Time (s)")
    plt.ylabel("Remaining nuclei")
    plt.title("Analytical vs Monte Carlo Radioactive Decay")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # Save the graph
    import os
    os.makedirs("results", exist_ok=True)
    plt.savefig("results/decay_comparison.png", dpi=300)

    print("Initial nuclei:", initial_nuclei)
    print("Decay constant:", decay_constant, "s^-1")
    print("Theoretical half-life:", np.log(2) / decay_constant, "s")
    print("Final analytical population:", analytical_population[-1])
    print("Final Monte Carlo population:", mc_population[-1])
    print("Graph saved to results/decay_comparison.png")

    plt.show()


if __name__ == "__main__":
    main()
