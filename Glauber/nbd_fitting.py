"""Fit a Glauber+negative-binomial model to RefMult3 data."""

import pathlib as path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize


def load_glauber_events(output_dir):
    """Load and combine the Glauber worker files."""
    files = sorted(path.Path(output_dir).glob("worker_*.npz"))
    if not files:
        raise FileNotFoundError(f"No worker files found in {output_dir}")

    npart = []
    ncoll = []
    b = []
    for file in files:
        with np.load(file) as data:
            npart.append(data["Npart"])
            ncoll.append(data["Ncoll"])
            b.append(data["b"])

    return np.concatenate(npart), np.concatenate(ncoll), np.concatenate(b)


def generate_nbd_multiplicity(npart, ncoll, alpha, mu, k, seed=12345):
    """Generate one RefMult3 value for each Glauber event."""
    nsource = (1 - alpha) * npart / 2 + alpha * ncoll
    shape = np.maximum(k * nsource, np.finfo(float).eps)
    probability = k / (k + mu)
    rng = np.random.default_rng(seed)
    return rng.negative_binomial(shape, probability).astype(np.int32)


def normalized_histogram(values, maximum):
    """Make a normalized integer histogram with overflow in the last bin."""
    values = np.clip(values.astype(int), 0, maximum)
    counts = np.bincount(values, minlength=maximum + 1)
    return counts / counts.sum()


def fit_nbd_parameters(
    refmult3,
    npart,
    ncoll,
    seed=12345,
    max_glauber_events=80_000,
    maxiter=120,
):
    """Fit alpha, mu, and k to the measured RefMult3 distribution."""
    refmult3 = np.asarray(refmult3)
    refmult3 = refmult3[np.isfinite(refmult3) & (refmult3 >= 0)].astype(int)
    if len(refmult3) == 0:
        raise ValueError("RefMult3 is empty")

    npart = np.asarray(npart)
    ncoll = np.asarray(ncoll)
    if npart.shape != ncoll.shape:
        raise ValueError("Npart and Ncoll must have the same shape")

    if len(npart) > max_glauber_events:
        rng = np.random.default_rng(seed)
        keep = rng.choice(len(npart), max_glauber_events, replace=False)
        npart = npart[keep]
        ncoll = ncoll[keep]

    minimum = int(refmult3.min())
    maximum = max(20, int(np.percentile(refmult3, 99.9)) + 10)
    data_histogram = normalized_histogram(refmult3, maximum)

    def objective(parameters):
        alpha, mu, k = parameters
        simulated = generate_nbd_multiplicity(
            npart,
            ncoll,
            alpha,
            mu,
            k,
            seed + 1,
        )
        simulated = simulated[simulated >= minimum]
        if len(simulated) == 0:
            return np.inf

        simulated_histogram = normalized_histogram(simulated, maximum)
        variance = (
            data_histogram / len(refmult3)
            + simulated_histogram / len(simulated)
        )
        enough_events = (
            (data_histogram * len(refmult3) >= 5)
            & (simulated_histogram * len(simulated) >= 5)
            & (variance > 0)
        )
        if np.sum(enough_events) <= 4:
            return np.inf

        chi_squared = np.sum(
            (data_histogram[enough_events] - simulated_histogram[enough_events]) ** 2
            / variance[enough_events]
        )
        return chi_squared / (np.sum(enough_events) - 4)

    return minimize(
        objective,
        [0.1, 1.5, 1.0],
        method="Powell",
        bounds=((0, 1), (0.05, 10), (0.05, 20)),
        options={"maxiter": maxiter, "xtol": 1e-3, "ftol": 1e-5},
    )


def plot_fitted_distribution(refmult3, fitted_refmult3, output_file):
    """Plot the measured and fitted RefMult3 distributions."""
    refmult3 = np.asarray(refmult3).astype(int)
    fitted_refmult3 = np.asarray(fitted_refmult3).astype(int)
    fitted_refmult3 = fitted_refmult3[fitted_refmult3 >= refmult3.min()]

    maximum = int(
        max(
            np.percentile(refmult3, 99.9),
            np.percentile(fitted_refmult3, 99.9),
        )
    )
    bins = np.arange(maximum + 2) - 0.5
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

    for axis, log_scale in zip(axes, (False, True)):
        axis.hist(
            refmult3,
            bins=bins,
            density=True,
            histtype="step",
            linewidth=1.8,
            label="UrQMD RefMult3",
        )
        axis.hist(
            fitted_refmult3,
            bins=bins,
            density=True,
            histtype="step",
            linewidth=1.8,
            label="Glauber+NBD",
        )
        axis.set(
            xlabel="RefMult3",
            ylabel="Probability",
            yscale="log" if log_scale else "linear",
        )
        axis.legend()

    fig.tight_layout()
    fig.savefig(output_file, dpi=180, bbox_inches="tight")
    plt.close(fig)
