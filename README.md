# Multi-class Kinetic Traffic Flow Model

## Description
This project provides a numerical implementation of the multi-class kinetic traffic flow model based on a non-local Prigogine-Herman framework, as derived by C. Mezquita-Nieto, P. Goatin, and A. Klar.

The research presents a discrete-velocity kinetic framework that models interactions between different vehicle classes through braking and relaxation terms, solved either by using path-conservative finite volume methods for the general velocity case or by using the transformation to conservative variables for the 2-velocity case.

The implementation consists of two primary numerical schemes for 1D heterogeneous traffic flow:
1. **Path-Conservative Finite Volume Scheme:** Used to evaluate non-conservative matrix jump interactions at cell interfaces.
2. **FORCE Scheme:** A robust finite volume method (integrating Lax-Friedrichs and Lax-Wendroff fluxes) for conservative 2-velocity modeling.

This model is intended for research in traffic engineering, kinetic theory, and hyperbolic conservation laws.

## Correction (October 2026)

In the version used for the preprint of the article (release [`v1.0-article`](https://github.com/carmen-mezquita/multi-class-kinetic-traffic-flow-model/releases/tag/v1.0-article)), the interaction matrix in `Ac_block.m` was set to zero whenever the total density satisfied $r \geq R_c$. This is not the model of the article, where the interaction term is used for all $r < 1$ and the maximal class densities $R_c$ only enter through the equilibrium functions. It affected the results of Fig. 6.5 (Test-case 3, initial data with mean total density above $R_2$). All other results are unaffected.

The current version uses the interaction matrix for all $r < 1$. At the endpoints of the integration path with $r = 1$ (jam states), the matrix is evaluated at its limit along the path.

## Repository contents
Two main functions to consider within the project:
1. `discvelo_path_conservative.m`: Main simulation engine for the non-local kinetic model when considering more than two velocities.
2. `discvelo_two_velocities.m`: Main simulation engine for the conservative 2-velocity model, using a FORCE numerical flux.

Additionally, there are helper functions available, including initialization routines, boundary condition handling, flux calculation and visualization tools for decay of solutions and density profiles. The documentation for the different functions is written as a header for each file.

## Visuals
The simulation engines generate diagnostic visualizations, including:
* **`I_t.png`**: Decay of density perturbations over time.
* **`total_mass.png`**: Tracking of total mass conservation.
* **Real-time plotting**: The scripts output real-time updates of density ($\rho$) and flux ($q$) profiles.

## Installation
The model requires **MATLAB** (R2020a or later, for `exportgraphics`). No additional toolboxes are needed.

1. Clone or download the repository.
2. Add the repository and its subfolders to the MATLAB path:
   ```matlab
   addpath(genpath('path/to/multi-class-kinetic-traffic-flow-model'));
   ```

## Usage
All parameters are set in `param.m`, which builds the configuration struct `par`. The test cases of the article are selected with `par.ic` (initial condition) and `par.bc` (boundary condition).

**Example execution:**

```matlab
% 1. Load the default parameters (all fields of par)
run('param.m');

% 2. Change some of them, e.g.
par.M     = 20;   % M+1 discrete velocities
par.t_end = 50;   % final time
par.ic    = 1;    % initial condition (test case)
par.rel   = 0;    % 1: with relaxation, 0: without
% (the grid size nx_inner and eps are changed directly in param.m, since other fields are computed from them)

% 3. Run the simulation
[I, total_mass, timeline, rho1, rho2, q1, q2] = discvelo_path_conservative(par);
% or, for the 2-velocity model (M = 1):
[I, total_mass, timeline, rho1, rho2, q1, q2] = discvelo_two_velocities(par);
```
