function f = boundary_conditions(nx,nbc,f,bc,V)
    %% BOUNDARY_CONDITIONS
    %       Extrapolates ghost cell elements based on boundary 
    %       configurations.
    %
    %   SYNTAX:
    %       f = boundary_conditions(nx, nbc, f, bc, V)
    %
    %   DESCRIPTION:
    %       Applies specialized spatial boundary conditions across the 
    %       ghost-cell regions of the velocity distribution array. Supports 
    %       absorbing free conditions, periodic updates, and restricted 
    %       Fan-Work inflows.
    %   INPUTS:
    %       nx    - Spatial domain size.
    %       nbc   - Number of ghost/boundary cells.
    %       f     - State distribution matrix.
    %       bc    - Boundary condition index.
    %       V     - Velocity vector.
    %
    %   OUTPUTS:
    %       f     - Distribution matrix with updated boundary conditions.
    %
    %   See also: DISCVELO_PATH_CONSERVATIVE.
    %
    %% CODE:

    switch bc
        case 1 % Free boundary condition profile (Zero Neumann replication)
            f(:,1:nbc) = f(:,nbc+1);
            f(:,nx-nbc+1:nx) = f(:,end-nbc);  

        case 2 % Periodic loop-back wrapping
            f(:,1:nbc) = f(:,nx-2*nbc+1:nx-nbc);
            f(:,nx-nbc+1:nx) = f(:,nbc+1:2*nbc);

        case 3 % Fan-Work benchmark conditions (Absorbing inflow, free open outflow)
            f(:,1:nbc) = 0;
            f(:,nx-nbc+1:nx) = f(:,end-nbc);  
    end   
end