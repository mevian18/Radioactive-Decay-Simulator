
"""Statistical study of Monte Carlo radioactive decay."""

import numpy as np
import matplotlib.pyplot as plt

from src.monte_carlo import simulate_decay


# 1. Define the physical parameters
initial_nuclei = 1000
decay_constant = 0.1  # per second
total_time = 30.0     # seconds

# 2. Define the simulation settings
number_of_simulations = 1000
time_step = 0.1       # seconds
random_seed = 42




# 3. Calculate the analytical prediction
survival_probability = np.exp(-decay_constant * total_time)
analytical_population = initial_nuclei * survival_probability

# 4. Run independent Monte Carlo simulations
final_populations = np.zeros(number_of_simulations, dtype=int)

rng = np.random.default_rng(random_seed)

for i in range(number_of_simulations):
    # Give each simulation a reproducible, independent random seed
    simulation_seed = int(rng.integers(0, 2**32))

    times, populations = simulate_decay(
        initial_nuclei=initial_nuclei,
        decay_constant=decay_constant,
        time_step=time_step,
        total_time=total_time,
        seed=simulation_seed,
    )

    # Save the population remaining at the final time
    final_populations[i] = populations[-1]

# 5. Calculate the statistical results
mean_population = np.mean(final_populations)
std_population = np.std(final_populations, ddof=1)

print("\nMONTE CARLO STATISTICAL STUDY")
print("----------------------------")
print(f"Number of simulations: {number_of_simulations}")
print(f"Analytical prediction: {analytical_population:.3f}")
print(f"Monte Carlo mean: {mean_population:.3f}")
print(f"Monte Carlo standard deviation: {std_population:.3f}")
print(
    f"Difference from analytical prediction: "
    f"{mean_population - analytical_population:.3f}"
)



# 6. Calculate theoretical statistical quantities
theoretical_std = np.sqrt(
    initial_nuclei
    * survival_probability
    * (1 - survival_probability)
)

standard_error = std_population / np.sqrt(number_of_simulations)

print(f"Theoretical standard deviation: {theoretical_std:.3f}")
print(f"Standard error of Monte Carlo mean: {standard_error:.3f}")

# 7. Plot the distribution of final populations
plt.figure(figsize=(10, 6))

plt.hist(
    final_populations,
    bins=25,
    edgecolor="black",
    alpha=0.75,
)

plt.axvline(
    analytical_population,
    linestyle="--",
    linewidth=2,
    label=f"Analytical prediction: {analytical_population:.2f}",
)

plt.axvline(
    mean_population,
    linestyle="-",
    linewidth=2,
    label=f"Monte Carlo mean: {mean_population:.2f}",
)

plt.xlabel("Final number of surviving nuclei")
plt.ylabel("Frequency")
plt.title("Distribution of Monte Carlo Radioactive Decay Outcomes")
plt.legend()
plt.grid(True, alpha=0.25)
plt.tight_layout()

plt.savefig("results/monte_carlo_distribution.png", dpi=300)
plt.show()

print("Histogram saved to results/monte_carlo_distribution.png")
