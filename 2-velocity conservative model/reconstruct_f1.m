function [f1,rho1] = reconstruct_f1(u1,u2,m,nx,nbc,bc)
    %% RECONSTRUCT_F1 
    %       Computes kinetic variables for class 1 from cumulative
    %       parameters.
    %
    %   SYNTAX:
    %       [f1, rho1] = reconstruct_f1(u1, u2, m, nx, nbc, bc)
    %
    %   DESCRIPTION:
    %       Maps cumulative system variables to class 1 kinetic 
    %       distribution.
    %
    %   INPUTS:
    %       u1, u2 - Cumulative variable states.
    %       m      - Number of velocity classes.
    %       nx     - Spatial dimension.
    %       nbc    - Boundary cell count.
    %       bc     - Boundary condition index.
    %
    %   OUTPUTS:
    %       f1     - Class 1 kinetic distribution.
    %       rho1   - Class 1 macroscopic density.
    %
    %   See also: RECONSTRUCT_F2, DISCVELO_TWO_VELOCITIES.
    %
    %% CODE:

    % u1(1,:) is N_0^1, u1(2,:) is N_1^1
    % u2(1,:) is N_0^2, u2(2,:) is N_0^2
    f1(2,:) = u2(1,:)./u1(2,:);
    f1(1,:) = u1(1,:) - f1(2,:);

    f1 = boundary_conditions(nx,nbc,f1,bc);
    rho1(:) = f1(1,:) + f1(2,:);
end