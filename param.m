%% param.m
%       File that contains all the parameters used in the different
%       calculations.

%% ----------- Main parameters: -----------

par.N = 2; % number of vehicle classes taken
par.M = 1; % number of discretization choices for the velocities, v_0,...,v_M 
           % 1, 10, 20.

par.rel = 0; % rel=0, we do not do the relaxation step; 
             % rel=1, we do

par.ic = 3; % ic=1: Riemann problem, no relaxation (test case 1)
            % ic=2: Riemann problem, relaxation (test case 2)
            % ic=3: periodic initial conditions (test case 3)
            
par.bc = 1; % bc=1: free bc, 
            % bc=2: periodic bc
            % bc=3: test case 2 (0 inflow, free outflow)

par.eps = 0.01; % relaxation term: 1, 0.01, 0.001
par.T_eps = par.eps*[1 1]; % relaxation times vector
par.V_eps = [0.6 1]; % maximal velocity vector for fundamental diagram
par.R_eps = [1 0.6]; % maximal density vector for fundamental diagram

%% ----------- Discretization settings: -----------

par.nx_inner = 1000; % number of cells in the domain -> road discretization 
                     % 1000, 1500, 2000
                     
par.nbc = 1; % number of cells at the boundaries
par.nx = par.nx_inner + 2*par.nbc;      % total number of cells
par.inner = par.nbc+1 : par.nx-par.nbc; % indices of the inner cells
%par.inter = par.nbc : par.nx-par.nbc; 
%par.inter2 = par.nbc+1 : par.nx; %= par.inter+1

%% ----------- Spatial settings: -----------

par.x_low = 0.0; % lower endpoint of the domain
par.x_up  = 1.0; % upper endpoint of the domain 
                 % 1.0, 50.0;

par.dx = (par.x_up - par.x_low) / par.nx_inner;
par.x = par.x_low+par.dx/2 : par.dx : par.x_up-par.dx/2;
par.x_long = par.x_low-par.dx/2 : par.dx : par.x_up+par.dx/2;

%% ----------- Time settings: -----------

par.t_begin = 0.0; % initial time
par.t_end = 0.3;   % final time 
                   % 0.3, 10, 15, 20, 33, 50, 80;