
"""Study how background radiation affects activity measurements."""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from src.detector_model import estimate_activity


# --------------------------------------------------
# 1. Define the physical system
# --------------------------------------------------

TRUE_ACTIVITY = 100.0       # Source activity in Bq
EFFICIENCY = 0.80            # Detection efficiency
BACKGROUND_RATES = [1, 5, 10, 20]  # Background counts per second
COUNTING_TIMES = [1, 5, 10, 30, 60]  # Seconds
N_TRIALS = 2000

rng = np.random.default_rng(42)


# --------------------------------------------------
# 2. Run the Monte Carlo experiment
# --------------------------------------------------

results = []

for background_rate in BACKGROUND_RATES:
    for duration in COUNTING_TIMES:

        activity_estimates = []
        reported_uncertainties = []

        for _ in range(N_TRIALS):

            # Counts measured with source and background
            gross_counts = rng.poisson(
                (EFFICIENCY * TRUE_ACTIVITY + background_rate)
                * duration
            )

            # Independent background-only measurement
            background_counts = rng.poisson(
                background_rate * duration
            )

            # Estimate source activity and its uncertainty
            estimate = estimate_activity(
                gross_counts=gross_counts,
                background_counts=background_counts,
                source_time=duration,
                background_time=duration,
                efficiency=EFFICIENCY,
            )

            activity_estimates.append(
                estimate["estimated_activity_bq"]
            )

            reported_uncertainties.append(
                estimate["activity_uncertainty_bq"]
            )

        # Convert results to numerical arrays
        activity_estimates = np.array(activity_estimates)
        reported_uncertainties = np.array(reported_uncertainties)

        # Summarize the simulated measurements
        mean_activity = np.mean(activity_estimates)
        bias = mean_activity - TRUE_ACTIVITY
        empirical_sd = np.std(activity_estimates, ddof=1)
        mean_uncertainty = np.mean(reported_uncertainties)

        # Theoretical SD for equal source and background times
        theoretical_sd = np.sqrt(
            (EFFICIENCY * TRUE_ACTIVITY + 2 * background_rate)
            / (EFFICIENCY**2 * duration)
        )

        results.append({
            "background_rate": background_rate,
            "duration": duration,
            "mean_activity": mean_activity,
            "bias": bias,
            "empirical_sd": empirical_sd,
            "mean_uncertainty": mean_uncertainty,
            "theoretical_sd": theoretical_sd,
        })


# --------------------------------------------------
# 3. Display the results
# --------------------------------------------------

print("\nBACKGROUND-RATE STUDY")
print("=" * 100)
print(f"True source activity: {TRUE_ACTIVITY:.1f} Bq")
print(f"Detection efficiency: {EFFICIENCY:.2f}")
print(f"Trials per condition: {N_TRIALS}")
print()

print(
    f"{'Background':>10} {'Time (s)':>9} {'Mean (Bq)':>12} "
    f"{'Bias (Bq)':>11} {'Empirical SD':>14} "
    f"{'Mean uncertainty':>18} {'Theory SD':>12}"
)
print("-" * 100)

for row in results:
    print(
        f"{row['background_rate']:>10.1f} "
        f"{row['duration']:>9.0f} "
        f"{row['mean_activity']:>12.3f} "
        f"{row['bias']:>11.3f} "
        f"{row['empirical_sd']:>14.3f} "
        f"{row['mean_uncertainty']:>18.3f} "
        f"{row['theoretical_sd']:>12.3f}"
    )


# --------------------------------------------------
# 4. Plot uncertainty against background rate
# --------------------------------------------------

fig, ax = plt.subplots(figsize=(9, 6))

for duration in COUNTING_TIMES:
    selected = [
        row for row in results
        if row["duration"] == duration
    ]

    ax.plot(
        [row["background_rate"] for row in selected],
        [row["empirical_sd"] for row in selected],
        marker="o",
        label=f"{duration} s",
    )

ax.set_title("Effect of Background Rate on Activity Uncertainty")
ax.set_xlabel("Background rate (counts/s)")
ax.set_ylabel("Empirical SD of activity estimates (Bq)")
ax.grid(True, alpha=0.3)
ax.legend(title="Counting time")
fig.tight_layout()

output_dir = Path("results")
output_dir.mkdir(exist_ok=True)

plot_path = output_dir / "background_rate_uncertainty.png"
fig.savefig(plot_path, dpi=300)
plt.show()

print(f"\nGraph saved to: {plot_path}")
