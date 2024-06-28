# Mehrpooya2023_RandomWalk1DNetwork
Computer code in D used for discrete model simulations and Julia code for continuum model simulations presented in the research article _Mechanical cell interactions on curved interfaces_ by Pascal R Buenzli, Shahak Kuba, Ryan J Murphy, and Matthew J Simpson (2024). Preprint available at: https://arxiv.org/abs/TBA.

See instructions in main.d for how to run and visualise the discrete model.

## Simulation for Figure 3, 4, and 12
Run Model 4 with N=4, m=4, Hookean force in model_inputs.d;

Plot with n=1 coil periods per spring in snapshot.py => snapshot.pdf (Figs 3, 12);

Plot trajectories.py => trajectories.pdf (Fig 4a).

## Simulation for Figure 5
Run Model 0 with N=4, m=1, Hookean force, either straight_springs or curved_springs in model_inputs.d;

Plot with n=4 coil periods per spring in snapshot.py => snapshot.pdf.

## Simulation for Figures 6 and 7
Run Model 2 or 2c with N=8, m=1 or m=8, eitherh straight_springs or curved_springs in model_inputs.d;

Plot m=1 snapshots with n=4 coil periods per spring in snapshot.py => snapshot.pdf

## Simulations for Figure 8 and 9
Run Model 0 with N=8, m=1,2,4,or 8, t_end=2, and either hookean_force, lineardiff_force, or porousdiff_force in model_inputs.d;

Run spring-diffusion-continuum.jl after selecting the corresponding force law.

Plot: profiles.py => profiles.pdf (Fig 8) (uses linux command line tools sed, awk, tail, head)

