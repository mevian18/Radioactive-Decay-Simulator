"""Detector counting and activity estimation."""

import numpy as np
from statistics import NormalDist

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
# ---------------------------------------------------------------------------
# Detection limits (Currie, 1968)
# ---------------------------------------------------------------------------

def _z(probability):
    """Standard-normal quantile, e.g. _z(0.95) = 1.645."""
    return NormalDist().inv_cdf(probability)


def null_net_count_std(background_rate, source_time, background_time):
    """Std. dev. of net counts N = G - (Ts/Tb)B when NO source is present."""
    return np.sqrt(
        background_rate * source_time * (1.0 + source_time / background_time)
    )


def critical_level_counts(background_rate, source_time, background_time, alpha=0.05):
    """Currie critical level L_C (net counts).

    Declare 'source detected' when net counts exceed L_C. The false-positive
    probability is alpha. Gaussian approximation: valid when the expected
    background counts in the source measurement are roughly 10 or more.
    """
    sigma0 = null_net_count_std(background_rate, source_time, background_time)
    return _z(1.0 - alpha) * sigma0


def detection_limit_counts(
    background_rate, source_time, background_time, alpha=0.05, beta=0.05
):
    """Currie detection limit L_D (net counts).

    The true net count that is detected with probability 1 - beta by a test
    with false-positive probability alpha. For alpha = beta this reduces to
    L_D = k**2 + 2 * L_C.
    """
    sigma0 = null_net_count_std(background_rate, source_time, background_time)
    k_a, k_b = _z(1.0 - alpha), _z(1.0 - beta)
    l_c = k_a * sigma0
    return l_c + k_b**2 / 2.0 + (k_b / 2.0) * np.sqrt(
        k_b**2 + 4.0 * (l_c + sigma0**2)
    )


def minimum_detectable_activity(
    efficiency, background_rate, source_time, background_time, alpha=0.05, beta=0.05
):
    """Minimum detectable activity (Bq) = L_D / (efficiency * source_time)."""
    l_d = detection_limit_counts(
        background_rate, source_time, background_time, alpha, beta
    )
    return l_d / (efficiency * source_time)