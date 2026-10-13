"""Detector counting and activity estimation."""

import numpy as np


def estimate_activity(
    gross_counts,
    background_counts,
    source_time,
    background_time,
    efficiency,
):
    """Estimate source activity and its 1-sigma counting uncertainty.

    Accepts scalars or NumPy arrays of counts (one entry per Monte Carlo
    trial). Times are in seconds, efficiency is the probability of
    registering an event per source decay. The estimate is deliberately
    NOT clipped at zero, so it stays unbiased; negative values can occur
    at low counts. The uncertainty is unreliable below about 10 counts
    (it is exactly 0 when both counts are 0).
    """
    gross = np.asarray(gross_counts, dtype=float)
    back = np.asarray(background_counts, dtype=float)
    params = np.asarray([source_time, background_time, efficiency], dtype=float)

    if not (np.isfinite(gross).all() and np.isfinite(back).all()
            and np.isfinite(params).all()):
        raise ValueError("Inputs must be finite.")
    if (gross < 0).any() or (back < 0).any():
        raise ValueError("Counts cannot be negative.")
    if (gross != np.round(gross)).any():
        raise ValueError("Gross counts must be integers.")
    if (back != np.round(back)).any():
        raise ValueError("Background counts must be integers.")
    if source_time <= 0 or background_time <= 0:
        raise ValueError("Measurement times must be positive.")
    if not 0 < efficiency <= 1:
        raise ValueError("Efficiency must be greater than 0 and at most 1.")

    gross_rate = gross / source_time
    background_rate = back / background_time
    net_rate = gross_rate - background_rate
    rate_uncertainty = np.sqrt(gross / source_time**2 + back / background_time**2)

    return {
        "gross_rate": gross_rate,
        "background_rate": background_rate,
        "net_rate": net_rate,
        "estimated_activity_bq": net_rate / efficiency,
        "activity_uncertainty_bq": rate_uncertainty / efficiency,
    }


def theoretical_activity_std(
    true_activity, efficiency, background_rate, source_time, background_time
):
    """Exact standard deviation of the activity estimator (Poisson theory)."""
    variance = (
        (efficiency * true_activity + background_rate) / source_time
        + background_rate / background_time
    ) / efficiency**2
    return np.sqrt(variance)