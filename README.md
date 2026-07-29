# Deep_Learning_Npart_Predictor

This repository develops a deep learning workflow for predicting the number of participant nucleons, $\(N_{\mathrm{part}}\)$, on an event-by-event basis in heavy-ion collision simulations.

The goal is to use particle-level event information from UrQMD to estimate $\(N_{\mathrm{part}}\)$ more directly than standard multiplicity-based centrality cuts. A better event-by-event estimate of $\(N_{\mathrm{part}}\)$ may help reduce volume fluctuations in multiplicity distributions, especially when studying observables such as proton cumulants.

## Glauber event generator

The [`glauber`](glauber) directory provides a configurable Python Monte Carlo
Glauber model that writes event-level `Npart`, `Ncoll`, multiplicity, impact
parameter, and participant geometry to numpy NPZ file.
