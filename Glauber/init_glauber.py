import nuclei as nc
import numpy as np

batch = 10000
cross_section_mb = 42.0
b_max = 15.0

nucleus1 = nc.Nucleus(A=197, R=6.38, a=0.535, batch_size=batch)
nucleus2 = nc.Nucleus(A=197, R=6.38, a=0.535, batch_size=batch)
Npart_trials, Ncoll_trials = nc.collision_geometry(
    cross_section_mb=cross_section_mb,
    b_max=b_max,
    nucleus1=nucleus1,
    nucleus2=nucleus2,
)

# Trials with no nucleon-nucleon collisions are not nucleus-nucleus events.
collision_mask = Ncoll_trials > 0
Npart = Npart_trials[collision_mask]
Ncoll = Ncoll_trials[collision_mask]

np.savez(
    "glauber_output.npz",
    Npart=Npart,
    Ncoll=Ncoll,
    x1=nucleus1.x[collision_mask],
    y1=nucleus1.y[collision_mask],
    x2=nucleus2.x[collision_mask],
    y2=nucleus2.y[collision_mask],
)
