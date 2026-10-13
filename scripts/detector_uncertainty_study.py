
"""Study the effect of counting time on activity uncertainty."""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from src.detector_model import estimate_activity


# 1. Define the simulated detector conditions
true_activity = 100.0       # Bq
efficiency = 0.8
background_rate = 5.0       # counts per second

# Measurement durations to investigate
counting_times = [1, 5, 10, 30, 60]  # seconds

# Number of independent experiments per duration
number_of_trials = 2000

# Reproducibility
rng = np.random.default_rng(42)

results = []

# 2. Run repeated source and background measurements
for duration in counting_times:

    activity_estimates = []
    estimated_uncertainties = []

    for _ in range(number_of_trials):

        # Source-plus-background counts
        gross_mean = (
            efficiency * true_activity + background_rate
        ) * duration

        gross_counts = rng.poisson(gross_mean)

        # Independent background measurement of equal duration
        background_mean = background_rate * duration
        background_counts = rng.poisson(background_mean)

        estimate = estimate_activity(
            gross_counts=gross_counts,
            background_counts=background_counts,
            source_time=duration,
            background_time=duration,
            efficiency=efficiency,
        )

        activity_estimates.append(
            estimate["estimated_activity_bq"]
        )
        estimated_uncertainties.append(
            estimate["activity_uncertainty_bq"]
        )

    # 3. Summarize results for this duration
    activity_estimates = np.asarray(activity_estimates)
    estimated_uncertainties = np.asarray(
        estimated_uncertainties
    )

    mean_estimate = np.mean(activity_estimates)
    bias = mean_estimate - true_activity
    empirical_std = np.std(activity_estimates, ddof=1)
    mean_reported_uncertainty = np.mean(
        estimated_uncertainties
    )
    relative_error_percent = (
        100 * abs(bias) / true_activity
    )

    results.append({
        "duration": duration,
        "mean": mean_estimate,
        "bias": bias,
        "std": empirical_std,
        "uncertainty": mean_reported_uncertainty,
        "relative_error": relative_error_percent,
    })


# 4. Print the results
print("COUNTING-TIME UNCERTAINTY STUDY")
print("=" * 78)
print(f"True source activity: {true_activity:.1f} Bq")
print(f"Detection efficiency: {efficiency:.2f}")
print(f"Background rate: {background_rate:.1f} counts/s")
print(f"Trials per duration: {number_of_trials}")
print()

print(
    f"{'Time (s)':>9} {'Mean (Bq)':>12} {'Bias (Bq)':>12} "
    f"{'Empirical SD':>14} {'Mean uncertainty':>18} "
    f"{'Abs. bias (%)':>14}"
)
print("-" * 84)

for row in results:
    print(
        f"{row['duration']:>9} "
        f"{row['mean']:>12.3f} "
        f"{row['bias']:>12.3f} "
        f"{row['std']:>14.3f} "
        f"{row['uncertainty']:>18.3f} "
        f"{row['relative_error']:>14.3f}"
    )


# 5. Plot empirical variation and theoretical uncertainty
durations = np.array([row["duration"] for row in results])
empirical_std = np.array([row["std"] for row in results])

theoretical_std = np.sqrt(
    (efficiency * true_activity + 2 * background_rate)
    / (efficiency**2 * durations)
)

plt.figure(figsize=(9, 6))

plt.plot(
    durations,
    empirical_std,
    "o-",
    label="Empirical SD across trials",
)

plt.plot(
    durations,
    theoretical_std,
    "s--",
    label="Theoretical SD",
)

plt.xlabel("Counting duration (s)")
plt.ylabel("Activity uncertainty (Bq)")
plt.title("Effect of Counting Duration on Activity Uncertainty")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

results_folder = Path("results")
results_folder.mkdir(exist_ok=True)

output_path = results_folder / "counting_time_uncertainty.png"
plt.savefig(output_path, dpi=300)
plt.show()

print()
print(f"Graph saved to: {output_path}")
