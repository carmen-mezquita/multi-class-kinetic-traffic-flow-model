function [f1,f2] = density_conservation(U_pc,f1_old,f2_old,V,lmbd,inner,m,nx)
    %% DENSITY_CONSERVATION 
    %       Enforces exact macroscopic mass conservation via projection 
    %       scaling.
    %
    %   SYNTAX:
    %       [f1, f2] = density_conservation(U_pc, f1_old, f2_old, V, lmbd, inner, m, nx)
    %
    %   DESCRIPTION:
    %       Projects non-conservative scheme steps back into exact macro 
    %       constraints by comparing updates against standard conservative 
    %       finite volume calculations.
    %
    %   INPUTS:
    %       U_pc                  - Post-path-conservative update state.
    %       f1_old, f2_old        - Pre-update distributions.
    %       V, lmbd, inner, m, nx - Simulation parameters.
    %
    %   OUTPUTS:
    %       f1, f2                - Density-constrained distributions.
    %
    %   See also: DISCVELO_PATH_CONSERVATIVE.
    %
    %% CODE:

    q1_old = V*f1_old; 
    q2_old = V*f2_old;

    rho1_old = sum(f1_old,1); 
    rho2_old = sum(f2_old,1);

    [rho1_target,rho2_target] = deal(zeros(1,nx));
    
    rho1_target(inner) = 0.5 * (rho1_old(inner-1) + rho1_old(inner+1)) ...
                         - (lmbd/2) * (q1_old(inner+1) - q1_old(inner-1));
    rho2_target(inner) = 0.5 * (rho2_old(inner-1) + rho2_old(inner+1)) ...
                         - (lmbd/2) * (q2_old(inner+1) - q2_old(inner-1));
    
    f1_pc = max(0,U_pc(1:m,:)); 
    f2_pc = max(f1_pc,U_pc(m+1:end,:));

    rho1_pc = sum(f1_pc,1); 
    rho2_pc = sum(f2_pc,1);
    
    f1 = f1_pc; 
    f2 = f2_pc; 
    
    f1(:,inner) = f1_pc(:,inner).*(rho1_target(inner)./max(rho1_pc(inner),1e-12));
    f2(:,inner) = f2_pc(:,inner).*(rho2_target(inner)./max(rho2_pc(inner),1e-12));
end