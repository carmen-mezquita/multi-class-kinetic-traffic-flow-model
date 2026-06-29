function [f2,rho2] = reconstruct_f2(u2,m,nx,nbc,bc)
    %% RECONSTRUCT_F2 
    %       Computes kinetic variables for class 2 from cumulative
    %       parameters.
    %
    %   SYNTAX:
    %       [f2, rho2] = reconstruct_f2(u2, m, nx, nbc, bc)
    %
    %   DESCRIPTION: 
    %       Maps cumulative system variables to class 2 kinetic 
    %       distribution.
    %
    %   INPUTS:
    %       u2     - Cumulative variable state for class 2.
    %       m      - Number of velocity classes.
    %       nx     - Spatial dimension.
    %       nbc    - Boundary cell count.
    %       bc     - Boundary condition index.
    %
    %   OUTPUTS:
    %       f2     - Class 2 kinetic distribution.
    %       rho2   - Class 2 macroscopic density.
    %
    %   See also: RECONSTRUCT_F1, DISCVELO_TWO_VELOCITIES.
    %
    %% CODE:

    % u1(1,:) is N_0^1, u1(2,:) is N_1^1
    % u2(1,:) is N_0^2, u2(2,:) is N_0^2
    f2(1,:) = 1 - (u2(1,:)./u2(2,:));
    f2(2,:) = u2(1,:).*(1./u2(2,:)-1);
    
    f2 = boundary_conditions(nx,nbc,f2,bc);
    rho2(:) = f2(1,:) + f2(2,:); 
end