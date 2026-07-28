import numpy as np
import pathlib as path
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import nbinom
"""
The nucleus class is just to contain the nucleon geometry in batches
with a shape (1000(default),x(Nucleon number),y(Nucleon number))
"""
class Nucleus:

    def __init__(self, A, R, a, batch_size=1000, b_max=15.0):
        self.A = A  # Mass number (protons + neutrons)
        self.R = R # woodsaxon radius parameter
        self.a = a # woodsaxon skin depth parameter
        self.b_max = b_max # maximum impact parameter
        self.batch_size = batch_size

    def generate_nucleon_positions(self):
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
        theta = np.arccos(1 - 2 * np.random.rand(n_samples)).astype(np.float32)
        phi = (2 * np.pi * np.random.rand(n_samples)).astype(np.float32)
        r = np.interp(
            np.random.random(n_samples),
            self.r_cdf(),
            np.linspace(0, self.b_max, 10_000),
        ).astype(np.float32)
        return r, theta, phi
    #Now just convert to cartesian coordinates
    def cartesian_coordinate_sample(self):
        r, theta, phi = self.polar_coordinate_sample(self.A * self.batch_size)
        x = (r * np.sin(theta) * np.cos(phi)).reshape(self.batch_size, self.A)
        y = (r * np.sin(theta) * np.sin(phi)).reshape(self.batch_size, self.A)
        return x, y

"""
The following methods compute Npart, Ncoll and b and write them to
a .npz file in the output folder
"""

def collision_geometry(
    cross_section_mb,
    b_max,
    nucleus1: Nucleus,
    nucleus2: Nucleus,
    event_chunk_size=64,
):
    if event_chunk_size <= 0:
        raise ValueError("event_chunk_size must be positive")

    nucleus1.generate_nucleon_positions()
    nucleus2.generate_nucleon_positions()

    batch_size = nucleus1.batch_size
    cross_section_fm2 = cross_section_mb / 10

    # Shift the nuclei based on the impact parameter
    b = (b_max * np.sqrt(np.random.random(batch_size))).astype(np.float32)
    x1 = nucleus1.x + b[:, None] / 2
    x2 = nucleus2.x - b[:, None] / 2

    Npart = np.empty(batch_size, dtype=np.int32)
    Ncoll = np.empty(batch_size, dtype=np.int32)
    interaction_distance_squared = np.float32(cross_section_fm2 / np.pi)

    # Process a small group of events at a time so the pairwise A x A
    # distance arrays do not cover the entire output batch at once.
    for start in range(0, batch_size, event_chunk_size):
        stop = min(start + event_chunk_size, batch_size)

        distance_squared = (
            x1[start:stop, :, None] - x2[start:stop, None, :]
        )
        np.square(distance_squared, out=distance_squared)

        y_distance_squared = (
            nucleus1.y[start:stop, :, None]
            - nucleus2.y[start:stop, None, :]
        )
        np.square(y_distance_squared, out=y_distance_squared)
        distance_squared += y_distance_squared
        del y_distance_squared

        interaction_mask = distance_squared < interaction_distance_squared
        participants1 = np.any(interaction_mask, axis=2)
        participants2 = np.any(interaction_mask, axis=1)

        Npart[start:stop] = (
            participants1.sum(axis=1) + participants2.sum(axis=1)
        )
        Ncoll[start:stop] = interaction_mask.sum(axis=(1, 2))

        del distance_squared, interaction_mask, participants1, participants2

    return Npart, Ncoll, b

def generate_events(
    cross_section_mb,
    b_max,
    nucleus1: Nucleus,
    nucleus2: Nucleus,
    n_events,
    output_dir=None,
):
    if nucleus1.batch_size != nucleus2.batch_size:
        raise ValueError("Both nuclei must use the same batch size")
    if n_events <= 0:
        raise ValueError("n_events must be positive")

    if output_dir is None:
        output_dir = path.Path(__file__).resolve().parent / "glauber_output"
    else:
        output_dir = path.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    batch_size = nucleus1.batch_size
    n_batches = (n_events + batch_size - 1) // batch_size

    for i in range(n_batches):
        Npart, Ncoll, b = collision_geometry(
            cross_section_mb,
            b_max,
            nucleus1,
            nucleus2,
        )

        # Truncate the final batch if necessary.
        n_this_batch = min(
            batch_size,
            n_events - i * batch_size,
        )

        Npart = Npart[:n_this_batch]
        Ncoll = Ncoll[:n_this_batch]
        b = b[:n_this_batch]

        # Remove trials that produced no collision.
        collision_mask = Ncoll > 0
        Npart = Npart[collision_mask]
        Ncoll = Ncoll[collision_mask]
        b = b[collision_mask]

        output_file = output_dir / f"batch_{i:05d}.npz"

        np.savez_compressed(
            output_file,
            Npart=Npart,
            Ncoll=Ncoll,
            b=b,
        )

    Npart_batches = []
    Ncoll_batches = []
    b_batches = []

    for i in range(n_batches):
        batch_file = output_dir / f"batch_{i:05d}.npz"

        with np.load(batch_file) as batch:
            Npart_batches.append(batch["Npart"])
            Ncoll_batches.append(batch["Ncoll"])
            b_batches.append(batch["b"])

    Npart = np.concatenate(Npart_batches)
    Ncoll = np.concatenate(Ncoll_batches)
    b = np.concatenate(b_batches)

    combined_output_file = output_dir / f"combined_{n_events}_events.npz"
    np.savez_compressed(
        combined_output_file,
        Npart=Npart,
        Ncoll=Ncoll,
        b=b,
    )

    for i in range(n_batches):
        batch_file = output_dir / f"batch_{i:05d}.npz"
        batch_file.unlink()

    return combined_output_file


"""
The following mehtods generate the glauber multiplicity and compute the fitted paramters
from a given multiplicity distribution.
-N_source is the number of NBD samples per event
"""

def get_Nsource(Npart, Ncoll, alpha=0.1):
    return (1 - alpha) * Npart / 2 + alpha * Ncoll

def generate_event_multiplicity(Nsource, k, mu):
    p = k / (k + mu)
    return np.random.negative_binomial(n=k*Nsource, p=p, size=Nsource.shape)

def mul(Npart, Ncoll, k, mu, alpha):
    Nsource = get_Nsource(Npart, Ncoll, alpha)
    return generate_event_multiplicity(Nsource, k, mu)

def chisquared(mult_hist, glauber_mul_hist):
    return np.sum((np.log(mult_hist) - np.log(glauber_mul_hist))**2)

def fit_parameters(Npart, Ncoll, mult_hist, initial_guess=(0.1, 1.0, 0.1)):
    def objective(params):
        k, mu, alpha = params
        glauber_mul_hist = mul(Npart, Ncoll, k, mu, alpha)
        return chisquared(mult_hist, glauber_mul_hist)

    result = minimize(objective, initial_guess, bounds=[(0.01, None), (0.01, None), (0.0, 1.0)])
    return result.x
