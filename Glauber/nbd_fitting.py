"""Fit a two-component Glauber+NBD model to an observed RefMult3 histogram.

This module intentionally contains all multiplicity-model code.  The companion
``nuclei.py`` module is limited to nuclear geometry and Glauber event generation.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize


@dataclass(frozen=True)
class NBDFitResult:
    alpha: float
    mu: float
    k: float
    objective: float
    success: bool


@dataclass(frozen=True)
class DistributionFitStatistics:
    """Goodness-of-fit diagnostics for normalized multiplicity histograms."""

    probability_mae: float
    count_mae: float
    ks_distance: float
    poisson_deviance_per_ndf: float
    pearson_chi2_per_ndf: float
    ndf: int
    fit_minimum: int
    fit_maximum: int
    observed_events: int
    simulated_events: int


def load_glauber_events(output_dir: str | Path) -> tuple[np.ndarray, ...]:
    """Load and concatenate all ``worker_*.npz`` Glauber outputs."""
    files = sorted(Path(output_dir).glob("worker_*.npz"))
    if not files:
        raise FileNotFoundError(f"No worker_*.npz files found in {output_dir}")

    npart_parts: list[np.ndarray] = []
    ncoll_parts: list[np.ndarray] = []
    b_parts: list[np.ndarray] = []
    for file in files:
        with np.load(file, allow_pickle=False) as data:
            npart_parts.append(data["Npart"])
            ncoll_parts.append(data["Ncoll"])
            b_parts.append(data["b"])
    return (
        np.concatenate(npart_parts),
        np.concatenate(ncoll_parts),
        np.concatenate(b_parts),
    )


def get_nsource(
    npart: np.ndarray,
    ncoll: np.ndarray,
    alpha: float,
) -> np.ndarray:
    """Two-component source count: (1-alpha) Npart/2 + alpha Ncoll."""
    return (1.0 - alpha) * np.asarray(npart) / 2.0 + alpha * np.asarray(ncoll)


def generate_nbd_multiplicity(
    npart: np.ndarray,
    ncoll: np.ndarray,
    alpha: float,
    mu: float,
    k: float,
    *,
    seed: int = 12345,
) -> np.ndarray:
    """Draw the total multiplicity from the source-summed NBD."""
    if not 0.0 <= alpha <= 1.0 or mu <= 0.0 or k <= 0.0:
        raise ValueError("Require 0 <= alpha <= 1 and positive mu and k")
    nsource = get_nsource(npart, ncoll, alpha)
    shape = np.maximum(k * nsource, np.finfo(float).eps)
    probability = k / (k + mu)
    rng = np.random.default_rng(seed)
    return rng.negative_binomial(shape, probability).astype(np.int32)


def normalized_histogram(values: np.ndarray, max_multiplicity: int) -> np.ndarray:
    """Return a normalized integer histogram with overflow in the final bin."""
    clipped = np.clip(np.asarray(values, dtype=np.int64), 0, max_multiplicity)
    counts = np.bincount(clipped, minlength=max_multiplicity + 1).astype(float)
    return counts / counts.sum()


def distribution_fit_statistics(
    observed_multiplicity: np.ndarray,
    simulated_multiplicity: np.ndarray,
    *,
    n_parameters: int = 3,
    minimum_expected_count: float = 5.0,
) -> DistributionFitStatistics:
    """Compare data and fitted multiplicity with standard binned diagnostics.

    The simulated sample is conditioned to the observed fit range before both
    histograms are normalized.  The final bin contains overflow.  Pearson
    chi-square and Poisson deviance use bins with at least
    ``minimum_expected_count`` expected entries.
    """
    observed = np.asarray(observed_multiplicity)
    simulated = np.asarray(simulated_multiplicity)
    observed = observed[np.isfinite(observed)].astype(np.int64)
    simulated = simulated[np.isfinite(simulated)].astype(np.int64)
    if observed.size == 0 or simulated.size == 0:
        raise ValueError("Observed and simulated multiplicity samples must be nonempty")

    fit_minimum = int(observed.min())
    fit_maximum = max(20, int(np.ceil(np.percentile(observed, 99.9))) + 10)
    observed = observed[observed >= fit_minimum]
    simulated = simulated[simulated >= fit_minimum]
    if simulated.size == 0:
        raise ValueError("No simulated events survive the observed fit range")

    observed_counts = np.bincount(
        np.clip(observed, 0, fit_maximum),
        minlength=fit_maximum + 1,
    ).astype(float)
    simulated_counts = np.bincount(
        np.clip(simulated, 0, fit_maximum),
        minlength=fit_maximum + 1,
    ).astype(float)
    observed_probability = observed_counts / observed_counts.sum()
    simulated_probability = simulated_counts / simulated_counts.sum()
    expected_counts = simulated_probability * observed_counts.sum()

    valid = expected_counts >= minimum_expected_count
    observed_valid = observed_counts[valid]
    expected_valid = expected_counts[valid]
    pearson = np.sum(
        (observed_valid - expected_valid) ** 2 / expected_valid
    )
    deviance_terms = expected_valid - observed_valid
    positive = observed_valid > 0
    deviance_terms[positive] = (
        observed_valid[positive]
        * np.log(observed_valid[positive] / expected_valid[positive])
        - observed_valid[positive]
        + expected_valid[positive]
    )
    deviance = 2.0 * np.sum(deviance_terms)
    # Normalized histograms have one constraint in addition to fitted params.
    ndf = max(int(valid.sum()) - n_parameters - 1, 1)

    return DistributionFitStatistics(
        probability_mae=float(
            np.mean(np.abs(observed_probability - simulated_probability))
        ),
        count_mae=float(np.mean(np.abs(observed_counts - expected_counts))),
        ks_distance=float(
            np.max(
                np.abs(
                    np.cumsum(observed_probability)
                    - np.cumsum(simulated_probability)
                )
            )
        ),
        poisson_deviance_per_ndf=float(deviance / ndf),
        pearson_chi2_per_ndf=float(pearson / ndf),
        ndf=ndf,
        fit_minimum=fit_minimum,
        fit_maximum=fit_maximum,
        observed_events=int(observed.size),
        simulated_events=int(simulated.size),
    )


def fit_nbd_parameters(
    observed_refmult3: np.ndarray,
    npart: np.ndarray,
    ncoll: np.ndarray,
    *,
    seed: int = 12345,
    max_glauber_events: int = 80_000,
    maxiter: int = 120,
) -> NBDFitResult:
    """Fit alpha, mu, and k with a two-sample Pearson histogram statistic."""
    observed = np.asarray(observed_refmult3)
    observed = observed[np.isfinite(observed) & (observed >= 0)].astype(np.int32)
    if observed.size == 0:
        raise ValueError("The observed RefMult3 sample is empty")

    npart = np.asarray(npart)
    ncoll = np.asarray(ncoll)
    if npart.shape != ncoll.shape:
        raise ValueError("Npart and Ncoll must have the same shape")

    rng = np.random.default_rng(seed)
    if len(npart) > max_glauber_events:
        keep = rng.choice(len(npart), size=max_glauber_events, replace=False)
        npart_fit = npart[keep]
        ncoll_fit = ncoll[keep]
    else:
        npart_fit = npart
        ncoll_fit = ncoll

    max_multiplicity = max(
        20,
        int(np.ceil(np.percentile(observed, 99.9))) + 10,
    )
    target_hist = normalized_histogram(observed, max_multiplicity)
    target_counts = target_hist * observed.size
    minimum_multiplicity = int(observed.min())

    def objective(parameters: np.ndarray) -> float:
        alpha, mu, k = parameters
        # Powell can evaluate slightly outside its bounds while bracketing.
        # Treat those trial points as invalid instead of aborting the fit.
        if not (0.0 <= alpha <= 1.0 and mu > 0.0 and k > 0.0):
            return float("inf")
        simulated = generate_nbd_multiplicity(
            npart_fit,
            ncoll_fit,
            alpha,
            mu,
            k,
            seed=seed + 1,
        )
        # Compare like with like: the notebook applies a RefMult3 trigger to
        # UrQMD, so the same lower threshold must condition the Glauber sample.
        simulated = simulated[simulated >= minimum_multiplicity]
        if simulated.size == 0:
            return float("inf")
        simulated_hist = normalized_histogram(simulated, max_multiplicity)
        simulated_counts = simulated_hist * simulated.size
        # Compare normalized histograms with the statistical variance from
        # both the UrQMD and finite Glauber Monte Carlo samples.
        variance = (
            target_hist / observed.size
            + simulated_hist / simulated.size
        )
        valid = (
            (target_counts >= 5.0)
            & (simulated_counts >= 5.0)
            & (variance > 0.0)
        )
        if valid.sum() <= 4:
            return float("inf")
        chi_square = np.sum(
            (target_hist[valid] - simulated_hist[valid]) ** 2
            / variance[valid]
        )
        ndf = valid.sum() - 4  # alpha, mu, k, and normalized-hist constraint
        return float(chi_square / ndf)

    starts = (
        np.array([0.10, 1.5, 1.0]),
        np.array([0.20, 1.0, 0.7]),
        np.array([0.05, 2.0, 2.0]),
    )
    bounds = ((0.0, 1.0), (0.05, 10.0), (0.05, 20.0))
    fits = [
        minimize(
            objective,
            start,
            method="Powell",
            bounds=bounds,
            options={"maxiter": maxiter, "xtol": 1e-3, "ftol": 1e-5},
        )
        for start in starts
    ]
    best = min(fits, key=lambda result: result.fun)
    alpha, mu, k = best.x
    return NBDFitResult(
        alpha=float(alpha),
        mu=float(mu),
        k=float(k),
        objective=float(best.fun),
        success=bool(best.success),
    )


def conditional_mean_calibration(
    multiplicity: np.ndarray,
    target: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Create E[target | integer multiplicity], interpolating empty bins."""
    multiplicity = np.asarray(multiplicity, dtype=np.int64)
    target = np.asarray(target, dtype=float)
    if multiplicity.shape != target.shape:
        raise ValueError("Multiplicity and target must have the same shape")
    max_multiplicity = int(multiplicity.max())
    counts = np.bincount(multiplicity, minlength=max_multiplicity + 1)
    sums = np.bincount(
        multiplicity,
        weights=target,
        minlength=max_multiplicity + 1,
    )
    grid = np.arange(max_multiplicity + 1, dtype=float)
    valid = counts > 0
    means = np.empty_like(grid)
    means[valid] = sums[valid] / counts[valid]
    means[~valid] = np.interp(grid[~valid], grid[valid], means[valid])
    return grid, means


def predict_from_calibration(
    observed_multiplicity: np.ndarray,
    grid: np.ndarray,
    conditional_mean: np.ndarray,
) -> np.ndarray:
    """Predict by interpolating a conditional-mean Glauber calibration."""
    return np.interp(
        np.asarray(observed_multiplicity, dtype=float),
        grid,
        conditional_mean,
        left=conditional_mean[0],
        right=conditional_mean[-1],
    )


def plot_fitted_distribution(
    observed_refmult3: np.ndarray,
    fitted_refmult3: np.ndarray,
    output_file: str | Path,
) -> None:
    """Save linear- and log-scale RefMult3 fit comparisons."""
    observed = np.asarray(observed_refmult3, dtype=np.int64)
    fitted = np.asarray(fitted_refmult3, dtype=np.int64)
    fitted = fitted[fitted >= observed.min()]
    if fitted.size == 0:
        raise ValueError("No fitted events survive the observed RefMult3 threshold")
    statistics = distribution_fit_statistics(observed, fitted)
    maximum = int(max(np.percentile(observed, 99.9), np.percentile(fitted, 99.9)))
    bins = np.arange(maximum + 2) - 0.5
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    for ax, log_y in zip(axes, (False, True)):
        ax.hist(
            observed,
            bins=bins,
            density=True,
            histtype="step",
            linewidth=1.8,
            label="URQMD RefMult3",
        )
        ax.hist(
            fitted,
            bins=bins,
            density=True,
            histtype="step",
            linewidth=1.8,
            label="Fitted Glauber+NBD",
        )
        ax.set(
            xlabel="RefMult3",
            ylabel="Probability",
            title=(
                "Glauber+NBD fit"
                + (" (log scale)" if log_y else "")
                + f"\nMAE={statistics.probability_mae:.3g}, "
                + f"KS={statistics.ks_distance:.3g}"
            ),
            yscale="log" if log_y else "linear",
        )
        ax.legend()
    fig.tight_layout()
    fig.savefig(output_file, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _load_refmult3(feature_archive: str | Path) -> np.ndarray:
    with np.load(feature_archive, allow_pickle=False) as data:
        feature_names = data["feature_names"].astype(str)
        matches = np.flatnonzero(feature_names == "refmult3")
        if len(matches) != 1:
            raise KeyError("Feature archive must contain exactly one refmult3 column")
        return data["X"][:, int(matches[0])]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glauber-dir", type=Path, required=True)
    parser.add_argument("--urqmd-features", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--fit-sample-size", type=int, default=80_000)
    parser.add_argument("--maxiter", type=int, default=120)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    observed = _load_refmult3(args.urqmd_features)
    observed = observed[observed > 0]
    npart, ncoll, b = load_glauber_events(args.glauber_dir)
    fit = fit_nbd_parameters(
        observed,
        npart,
        ncoll,
        seed=args.seed,
        max_glauber_events=args.fit_sample_size,
        maxiter=args.maxiter,
    )
    fitted = generate_nbd_multiplicity(
        npart,
        ncoll,
        fit.alpha,
        fit.mu,
        fit.k,
        seed=args.seed + 2,
    )
    np.savez_compressed(
        args.output_dir / "nbd_fit_results.npz",
        alpha=fit.alpha,
        mu=fit.mu,
        k=fit.k,
        objective=fit.objective,
        success=fit.success,
        fitted_refmult3=fitted,
        Npart=npart,
        Ncoll=ncoll,
        b=b,
    )
    plot_fitted_distribution(
        observed,
        fitted,
        args.output_dir / "nbd_fit_refmult3.png",
    )
    print(
        f"alpha={fit.alpha:.5f}, mu={fit.mu:.5f}, k={fit.k:.5f}, "
        f"objective={fit.objective:.6g}, success={fit.success}"
    )


if __name__ == "__main__":
    main()