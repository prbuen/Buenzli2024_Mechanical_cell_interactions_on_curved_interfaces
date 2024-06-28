# Mechanical cell interactions on curved interfaces
Computer code in D used for discrete model simulations and Julia code for continuum model simulations presented in the research article _Mechanical cell interactions on curved interfaces_ by Pascal R Buenzli, Shahak Kuba, Ryan J Murphy, and Matthew J Simpson (2024). Preprint available at: [https://arxiv.org/abs/2406.19197](https://arxiv.org/abs/2406.19197).

See instructions in src/main.d for how to run and visualise the discrete model.

## Simulations for Figures 3, 4, and 12
Run Model 4 with N=4, m=4, Hookean force in model_inputs.d;

Plot with n=1 coil periods per spring in snapshot.py => snapshot.pdf (Figs 3);
![snapshot](https://github.com/prbuen/Buenzli2024_Mechanical_cell_interactions_on_curved_interfaces/assets/54585460/68581819-59cc-4173-b875-5d30d5133e41)


Plot trajectories.py => trajectories.pdf (Fig 4a).
![trajectories](https://github.com/prbuen/Buenzli2024_Mechanical_cell_interactions_on_curved_interfaces/assets/54585460/8d344549-1e1f-4876-9b93-2c6bebe41060)


## Simulation for Figure 5
Run Model 0 with N=4, m=1, Hookean force, either straight_springs or curved_springs in model_inputs.d;

Plot with n=4 coil periods per spring in snapshot.py => snapshot.pdf.

## Simulation for Figures 6 and 7
Run Model 2 with N=8, m=1 or m=8, or Model 2c with N=7, m=1 or m=8, either straight_springs or curved_springs in model_inputs.d;

Plot m=1 snapshots with n=4 coil periods per spring in snapshot.py => snapshot.pdf

## Simulations for Figure 8 and 9
Run Model 0 with N=8, m=1,2,4,or 8, t_end=2, and either hookean_force, lineardiff_force, or porousdiff_force in model_inputs.d;

Run spring-diffusion-continuum.jl after selecting the corresponding force law.

Plot: profiles.py => profiles.pdf (Fig 8) (uses linux command line tools sed, awk, tail, head)
![profiles](https://github.com/prbuen/Buenzli2024_Mechanical_cell_interactions_on_curved_interfaces/assets/54585460/4d480cb5-d7c6-4e26-8a97-33bf1e804b7b)


