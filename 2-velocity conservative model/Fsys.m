function func = Fsys(u,V,m,nx) 
    %% FSYS  
    %       Flux system function
    %
    %   SYNTAX:
    %       func = Fsys(u, V, m, nx)
    %
    %   DESCRIPTION:
    %       Computes the physical flux function at each grid point and for
    %       each velocity.
    %
    %   INPUTS:
    %       u      - Cumulative state variables.
    %       V      - Velocity vector.
    %       m, nx  - Discretization parameters.
    %
    %   OUTPUTS:
    %       func   - Flux vector components.
    %
    % See also: LAMBDA
    %
    %% CODE:

    u1 = u(1:2,:);
    u2 = u(3:4,:);

    func(1,:) = u2(1,:)./u1(2,:);
    func(2,:) = u1(2,:);
    func(3:4,:) = lambda(u2,V,m,nx).*u2; % u = N_i, func = lambda_i*N_i
end