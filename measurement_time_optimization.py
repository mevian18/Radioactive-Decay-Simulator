"""Optimize time allocation for radioactive activity measurements."""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from src.detector_model import estimate_activity, theoretical_activity_std

TRUE_ACTIVITY = 100.0
EFFICIENCY = 0.80
BACKGROUND_RATE = 5.0
TOTAL_TIME = 60.0
N_TRIALS = 20_000

rng = np.random.default_rng(42)

gross_rate = EFFICIENCY * TRUE_ACTIVITY + BACKGROUND_RATE
optimal_fraction = np.sqrt(gross_rate) / (np.sqrt(gross_rate) + np.sqrt(BACKGROUND_RATE))

# Fine grid that also contains the theoretical optimum itself
fractions = np.sort(np.append(np.linspace(0.05, 0.95, 19), optimal_fraction))


def theory_sd(fraction):
    ts = TOTAL_TIME * fraction
    return theoretical_activity_std(
        TRUE_ACTIVITY, EFFICIENCY, BACKGROUND_RATE, ts, TOTAL_TIME - ts
    )


rows = []
for f in fractions:
    ts, tb = TOTAL_TIME * f, TOTAL_TIME * (1 - f)
    gross = rng.poisson(gross_rate * ts, N_TRIALS)
    back = rng.poisson(BACKGROUND_RATE * tb, N_TRIALS)
    est = estimate_activity(gross, back, ts, tb, EFFICIENCY)["estimated_activity_bq"]
    sd = est.std(ddof=1)
    rows.append((f, ts, tb, est.mean(), sd, sd / np.sqrt(2 * (N_TRIALS - 1)), theory_sd(f)))

print(f"Theoretical optimum: source {optimal_fraction:.3f}, background {1-optimal_fraction:.3f}"
      f"  (Tb/Ts = {np.sqrt(BACKGROUND_RATE/gross_rate):.3f})")
print(f"{'Ts':>6} {'Tb':>6} {'Mean':>8} {'MC SD':>8} {'+/-':>6} {'Theory':>8}")
for f, ts, tb, m, sd, err, th in rows:
    print(f"{ts:6.1f} {tb:6.1f} {m:8.3f} {sd:8.3f} {err:6.3f} {th:8.3f}")

best = min(rows, key=lambda r: r[4])
print(f"\nMC minimum at source fraction {best[0]:.3f}")
print(f"SD at optimum {theory_sd(optimal_fraction):.3f} vs equal split {theory_sd(0.5):.3f} "
      f"-> {100*(1-theory_sd(optimal_fraction)/theory_sd(0.5)):.1f}% smaller")

x = [r[1] for r in rows]
fig, ax = plt.subplots(figsize=(9, 6))
ax.errorbar(x, [r[4] for r in rows], yerr=[r[5] for r in rows], fmt="o",
            capsize=3, label="Monte Carlo SD (+/- 1 s.e.)")
fine = np.linspace(0.05, 0.95, 400)
ax.plot(TOTAL_TIME * fine, theory_sd(fine), "--", label="Theoretical SD")
ax.axvline(TOTAL_TIME * optimal_fraction, linestyle=":",
           label=f"Theoretical optimum: {TOTAL_TIME*optimal_fraction:.1f} s")
ax.set_xlabel(f"Source measurement time (s)  [background time = {TOTAL_TIME:.0f} s minus this]")
ax.set_ylabel("Activity uncertainty (Bq)")
ax.set_title("Effect of Measurement-Time Allocation on Activity Uncertainty")
ax.grid(True, alpha=0.3)
ax.legend()
fig.tight_layout()

Path("results").mkdir(exist_ok=True)
plot_path = Path("results") / "measurement_time_optimization.png"
fig.savefig(plot_path, dpi=300)
plt.show()

print(f"\nGraph saved to: {plot_path}")