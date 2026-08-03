"""Nuclear geometry and Glauber event generation.

Multiplicity and NBD fitting live in ``nbd_fitting.py``.
"""

import pathlib as path

import numpy as np


class Nucleus:
    """Contain batched nucleon geometry for one nucleus."""

    def __init__(self, A, R, a, batch_size=1000, b_max=15.0):
        self.A = A
        self.R = R
        self.a = a
        self.b_max = b_max
        self.batch_size = batch_size

    def generate_nucleon_positions(self):
        self.x, self.y = self.cartesian_coordinate_sample()

    def woodsaxon_density(self, r, R, a):
        return 1 / (1 + np.exp((r - R) / a))

    def r_cdf(self):
        r_grid = np.linspace(0, self.b_max, 10_000)
        rho = self.woodsaxon_density(r_grid, self.R, self.a)
        pdf = 4 * np.pi * r_grid**2 * rho
        return np.cumsum(pdf / np.sum(pdf))

    def polar_coordinate_sample(self, n_samples):
        theta = np.arccos(1 - 2 * np.random.rand(n_samples)).astype(np.float32)
        phi = (2 * np.pi * np.random.rand(n_samples)).astype(np.float32)
        r = np.interp(
            np.random.random(n_samples),
            self.r_cdf(),
            np.linspace(0, self.b_max, 10_000),
        ).astype(np.float32)
        return r, theta, phi

    def cartesian_coordinate_sample(self):
        r, theta, phi = self.polar_coordinate_sample(self.A * self.batch_size)
        x = (r * np.sin(theta) * np.cos(phi)).reshape(self.batch_size, self.A)
        y = (r * np.sin(theta) * np.sin(phi)).reshape(self.batch_size, self.A)
        return x, y


def collision_geometry(
    cross_section_mb,
    b_max,
    nucleus1: Nucleus,
    nucleus2: Nucleus,
    event_chunk_size=64,
):
    """Compute Npart, Ncoll, and b for one batch."""
    if event_chunk_size <= 0:
        raise ValueError("event_chunk_size must be positive")

    nucleus1.generate_nucleon_positions()
    nucleus2.generate_nucleon_positions()

    batch_size = nucleus1.batch_size
    cross_section_fm2 = cross_section_mb / 10
    b = (b_max * np.sqrt(np.random.random(batch_size))).astype(np.float32)
    x1 = nucleus1.x + b[:, None] / 2
    x2 = nucleus2.x - b[:, None] / 2

    npart = np.empty(batch_size, dtype=np.int32)
    ncoll = np.empty(batch_size, dtype=np.int32)
    interaction_distance_squared = np.float32(cross_section_fm2 / np.pi)

    for start in range(0, batch_size, event_chunk_size):
        stop = min(start + event_chunk_size, batch_size)
        distance_squared = x1[start:stop, :, None] - x2[start:stop, None, :]
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
        npart[start:stop] = participants1.sum(axis=1) + participants2.sum(axis=1)
        ncoll[start:stop] = interaction_mask.sum(axis=(1, 2))
        del distance_squared, interaction_mask, participants1, participants2

    return npart, ncoll, b


def generate_events(
    cross_section_mb,
    b_max,
    nucleus1: Nucleus,
    nucleus2: Nucleus,
    n_events,
    output_dir=None,
):
    """Generate Glauber geometry events and write one combined NPZ archive."""
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

    for index in range(n_batches):
        npart, ncoll, b = collision_geometry(
            cross_section_mb,
            b_max,
            nucleus1,
            nucleus2,
        )
        n_this_batch = min(batch_size, n_events - index * batch_size)
        npart = npart[:n_this_batch]
        ncoll = ncoll[:n_this_batch]
        b = b[:n_this_batch]

        collision_mask = ncoll > 0
        output_file = output_dir / f"batch_{index:05d}.npz"
        np.savez_compressed(
            output_file,
            Npart=npart[collision_mask],
            Ncoll=ncoll[collision_mask],
            b=b[collision_mask],
        )

    npart_batches = []
    ncoll_batches = []
    b_batches = []
    for index in range(n_batches):
        batch_file = output_dir / f"batch_{index:05d}.npz"
        with np.load(batch_file) as batch:
            npart_batches.append(batch["Npart"])
            ncoll_batches.append(batch["Ncoll"])
            b_batches.append(batch["b"])

    npart = np.concatenate(npart_batches)
    ncoll = np.concatenate(ncoll_batches)
    b = np.concatenate(b_batches)
    combined_output_file = output_dir / f"combined_{n_events}_events.npz"
    np.savez_compressed(
        combined_output_file,
        Npart=npart,
        Ncoll=ncoll,
        b=b,
    )

    for index in range(n_batches):
        (output_dir / f"batch_{index:05d}.npz").unlink()

    return combined_output_file