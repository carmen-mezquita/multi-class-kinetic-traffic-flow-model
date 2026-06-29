# Multi-class Kinetic Traffic Flow Model

## Name
Multi-class Kinetic Traffic Flow Model

## Description
This project provides a numerical implementation of the multi-class kinetic traffic flow model based on a non-local Prigogine-Herman framework, as derived by C. Mezquita-Nieto, P. Goatin, and A. Klar.

The research presents a discrete-velocity kinetic framework that models interactions between different vehicle classes through braking and relaxation terms, solved either by using path-conservative finite volume methods for the general velocity case or by using the transformation to conservative variables for the 2-velocity case.

The implementation consists of two primary numerical schemes for 1D heterogeneous traffic flow:
1. **Path-Conservative Finite Volume Scheme:** Used to evaluate non-conservative matrix jump interactions at cell interfaces.
2. **FORCE Scheme:** A robust finite volume method (integrating Lax-Friedrichs and Lax-Wendroff fluxes) for conservative 2-velocity modeling.

This model is intended for research in traffic engineering, kinetic theory, and hyperbolic conservation laws.

## Repository contents
Two main functions to consider within the project:
1. `discvelo_path_conservative.m`: Main simulation engine for the non-local kinetic model when considering more than two velocities.
2. `discvelo_two_velocities`: Main simulation engine for the conservative 2-velocity model, using a FORCE numerical flux.

Additionally, there are helper functions available, including initialization routines, boundary condition handling, flux calculation and visualization tools for decay of solutions and density profiles. The documentation for the different functions is written as a header for each file.

## Visuals
The simulation engines generate diagnostic visualizations, including:
* **`I_t.png`**: Decay of density perturbations over time.
* **`total_mass.png`**: Tracking of total mass conservation.
* **Real-time plotting**: The scripts output real-time updates of density ($\rho$) and flux ($q$) profiles.

## Installation
The model requires **MATLAB**, and the Optimization/Statistics Toolbox for the plotting and data saving functions. No complex installation or external toolboxes are strictly required for the core calculations, though standard plotting functions are utilized.

1. Clone the repository to your local machine.
2. Add the directory to your MATLAB path.
3. Ensure your current folder contains all supporting files (`boundary_conditions.m`, `initial_condition.m`, `flux_FORCE.m`, etc.).

## Usage
The main simulation functions are designed to be called with a configuration struct (`par`) that defines the spatial domain, discretization, and physical parameters (such as relaxation time or vehicle class capacity).

**Example execution for the path-conservative model:**

```matlab
% 1. Define simulation parameters
par.nx = 200;       % Spatial grid points
par.M = 1;          % Velocity classes
par.t_begin = 0;    % Start time
par.t_end = 50;     % End time
par.ic = 1;         % Test-case identifier

% 2. Run the simulation
[I, total_mass, timeline, rho1, rho2, q1, q2] = discvelo_path_conservative(par);

