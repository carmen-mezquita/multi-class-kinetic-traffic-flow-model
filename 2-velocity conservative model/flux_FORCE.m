function [flux1,flux2] = flux_FORCE(u1,v1,u2,v2,lmbd,V,m,nx)
    %% FLUX_FORCE  
    %       Computes the FORCE numerical flux (First-Order Centered Scheme 
    %       with Evolution).
    %
    %   SYNTAX:
    %       [flux1, flux2] = flux_FORCE(u1, v1, u2, v2, lmbd, V, m, nx)
    %
    %   DESCRIPTION: 
    %       Combines Lax-Friedrichs and Lax-Wendroff fluxes for robustness.
    %
    %   INPUTS:
    %       u1, u2       - Left state variables.
    %       v1, v2       - Right state variables.
    %       lmbd         - CFL coefficient (dt/dx).
    %       V            - Velocity vector.
    %       m, nx        - Discretization parameters.
    %
    %   OUTPUTS:
    %       flux1, flux2 - Interface fluxes.
    %
    % See also: FSYS, LAMBDA

    %% CODE:
    
    u = [u1; u2];
    v = [v1; v2];

    F_u(:,:) = Fsys(u(:,:),V,m,nx);
    F_v(:,:) = Fsys(v(:,:),V,m,nx);

    % Lax-Friedrichs flux:
    F_LF(:,:) = 0.5*(F_u(:,:) + F_v(:,:)) - (1/lmbd)*0.5*(v(:,:) - u(:,:));

    % Lax-Wendroff midpoint predictor: second-order accurate predictor of 
    % the solution at the midpoint in time and space.
    ustar(:,:) = (u(:,:) + v(:,:))/2 - (lmbd/2)*(F_v(:,:) - F_u(:,:));

    % Lax-Wendroff flux:
    F_LW(:,:) = Fsys(ustar(:,:),V,m,nx);

    % FORCE flux, written in a reorganized form:
    F_FORCE(:,:) = 0.5*(F_LF(:,:) + F_LW(:,:));
    flux1 = F_FORCE(1:2,:);
    flux2 = F_FORCE(3:4,:);
end