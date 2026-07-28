import nuclei as nc
import numpy as np

n_events = 10000
cross_section_mb = 28.0
b_max = 15.0

nucleus1 = nc.Nucleus(A=197, R=6.38, a=0.535)
nucleus2 = nc.Nucleus(A=197, R=6.38, a=0.535)

nc.generate_events(cross_section_mb,b_max,nucleus1,nucleus2,n_events)
