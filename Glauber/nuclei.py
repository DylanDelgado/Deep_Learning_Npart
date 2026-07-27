import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import nbinom

class Nucleus:

    def __init__(self, A, R, a, batch_size, b_max=15.0):
        self.A = A  # Mass number (protons + neutrons)
        self.R = R # woodsaxon radius parameter
        self.a = a # woodsaxon skin depth parameter
        self.batch_size = batch_size # number of events to simulate
        self.b_max = b_max # maximum impact parameter
        self.x, self.y = self.cartesian_coordinate_sample()

    # The following methods are used give xyz coordinates of nucleons in a nucleus
    #First create a function for the woodsaxon density.
    def woodsaxon_density(self, r, R, a):
        return 1 / (1 + np.exp((r - R) / a))
    # To sample r from this density we compute what the cdf is.
    def r_cdf(self):
        r_grid = np.linspace(0, self.b_max, 10_000)
        rho = self.woodsaxon_density(r_grid, self.R, self.a)
        pdf = 4 * np.pi * r_grid**2 * rho
        return np.cumsum(pdf/np.sum(pdf))
    #Sample polar coordinates (note to have uniform solid angle we sample costheta uniformly)
    def polar_coordinate_sample(self, n_samples):
        theta = np.arccos(1 - 2 * np.random.rand(n_samples))
        phi = 2 * np.pi * np.random.rand(n_samples)
        r = np.interp(np.random.random(n_samples), self.r_cdf(), np.linspace(0, self.b_max, 10_000))
        return r, theta, phi
    #Now just convert to cartesian coordinates
    def cartesian_coordinate_sample(self):
        r, theta, phi = self.polar_coordinate_sample(self.A * self.batch_size)
        x = (r * np.sin(theta) * np.cos(phi)).reshape(self.batch_size, self.A)
        y = (r * np.sin(theta) * np.sin(phi)).reshape(self.batch_size, self.A)
        return x, y

def collision_geometry(cross_section_mb, b_max, nucleus1: Nucleus, nucleus2: Nucleus):
    if nucleus1.batch_size != nucleus2.batch_size:
        raise ValueError("Both nuclei must have the same batch size.")
    if cross_section_mb <= 0:
        raise ValueError("cross_section_mb must be positive.")
    if b_max < 0:
        raise ValueError("b_max cannot be negative.")

    batch_size = nucleus1.batch_size
    cross_section_fm2 = cross_section_mb / 10

    # Shift the nuclei based on the impact parameter
    b = b_max * np.sqrt(np.random.random(batch_size)) 
    x1 = nucleus1.x + b[:, None] / 2
    x2 = nucleus2.x - b[:, None] / 2
    
    # Calculate distances between nucleons in the two nuclei
    distances = np.sqrt((x1[:,:,None] - x2[:,None,:])**2 +
                        (nucleus1.y[:,:,None] - nucleus2.y[:,None,:])**2) 

    # Determine which nucleons are within the interaction range
    interaction_mask = distances < np.sqrt(cross_section_fm2 / np.pi)
    
    # Count participants and collisions separately for each event
    participants1 = np.any(interaction_mask, axis=2)
    participants2 = np.any(interaction_mask, axis=1)
    Npart = participants1.sum(axis=1) + participants2.sum(axis=1)
    Ncoll = interaction_mask.sum(axis=(1, 2))
    
    return Npart, Ncoll

def get_Nsource(Npart, Ncoll, alpha=0.1):
    return (1 - alpha) * Npart / 2 + alpha * Ncoll

def generate_multiplicity(Nsource, k, mu):
    p = k / (k + mu)
    return np.random.negative_binomial(n=k*Nsource, p=p, size=Nsource.shape)


def fit_multiplicity(
    urqmd_multiplicity,
    Npart,
    Ncoll,
    initial_guess=None,
    maxiter=500,
):
    urqmd_multiplicity = np.asarray(urqmd_multiplicity)
    Npart = np.asarray(Npart, dtype=float)
    Ncoll = np.asarray(Ncoll, dtype=float)

    if urqmd_multiplicity.ndim != 1 or urqmd_multiplicity.size == 0:
        raise ValueError("urqmd_multiplicity must be a non-empty 1D array.")
    if (
        not np.all(np.isfinite(urqmd_multiplicity))
        or np.any(urqmd_multiplicity < 0)
        or np.any(urqmd_multiplicity != np.floor(urqmd_multiplicity))
    ):
        raise ValueError(
            "urqmd_multiplicity must contain finite, non-negative integers."
        )
    if Npart.ndim != 1 or Ncoll.ndim != 1 or Npart.size != Ncoll.size:
        raise ValueError("Npart and Ncoll must be 1D arrays of equal length.")
    if Npart.size == 0:
        raise ValueError("Npart and Ncoll cannot be empty.")
    if (
        not np.all(np.isfinite(Npart))
        or not np.all(np.isfinite(Ncoll))
        or np.any(Npart <= 0)
        or np.any(Ncoll <= 0)
    ):
        raise ValueError("Npart and Ncoll must contain finite, positive values.")
    if not isinstance(maxiter, (int, np.integer)) or maxiter <= 0:
        raise ValueError("maxiter must be a positive integer.")

    multiplicity_values, observed_counts = np.unique(
        urqmd_multiplicity.astype(np.int64),
        return_counts=True,
    )

    if initial_guess is None:
        initial_alpha = 0.1
        initial_Nsource = get_Nsource(Npart, Ncoll, initial_alpha)
        initial_mu = max(
            float(np.mean(urqmd_multiplicity) / np.mean(initial_Nsource)),
            1e-3,
        )
        initial_guess = (1.0, initial_mu, initial_alpha)
    else:
        initial_guess = np.asarray(initial_guess, dtype=float)
        if initial_guess.shape != (3,) or not np.all(np.isfinite(initial_guess)):
            raise ValueError("initial_guess must contain finite (k, mu, alpha).")
        if (
            initial_guess[0] <= 0
            or initial_guess[1] <= 0
            or not 0 <= initial_guess[2] <= 1
        ):
            raise ValueError(
                "initial_guess requires k > 0, mu > 0, and 0 <= alpha <= 1."
            )

    def negative_log_likelihood(parameters):
        k, mu, alpha = parameters
        Nsource = get_Nsource(Npart, Ncoll, alpha)
        total_k = k * Nsource
        p = k / (k + mu)

        component_logpmf = nbinom.logpmf(
            multiplicity_values[:, None],
            total_k[None, :],
            p,
        )
        mixture_logpmf = (
            logsumexp(component_logpmf, axis=1) - np.log(Npart.size)
        )
        return -np.dot(observed_counts, mixture_logpmf)

    result = minimize(
        negative_log_likelihood,
        x0=initial_guess,
        method="L-BFGS-B",
        bounds=((1e-8, None), (1e-8, None), (0.0, 1.0)),
        options={"maxiter": int(maxiter)},
    )

    result.k, result.mu, result.alpha = result.x
    return result
