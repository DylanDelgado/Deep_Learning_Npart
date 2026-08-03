import nuclei as nc
import numpy as np
import pathlib as path
import sys

n_events = 2000
cross_section_mb = 42.0
b_max = 18.0
worker_id = int(sys.argv[1])

np.random.seed(1000 + worker_id)

nucleus1 = nc.Nucleus(A=197, R=6.38, a=0.535)
nucleus2 = nc.Nucleus(A=197, R=6.38, a=0.535)

output_dir = (
    path.Path(__file__).resolve().parent
    / "glauber_output"
    / f"worker_{worker_id:02d}_temporary"
)

combined_file = nc.generate_events(
    cross_section_mb,
    b_max,
    nucleus1,
    nucleus2,
    n_events,
    output_dir,
)

worker_file = output_dir.parent / f"worker_{worker_id:02d}.npz"
combined_file.replace(worker_file)
output_dir.rmdir()