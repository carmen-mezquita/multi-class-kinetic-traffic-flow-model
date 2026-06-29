function U = path_conservative_update(U_old,I_pi,lmbd,bc,nx,nbc,inner)
    %% PATH_CONSERVATIVE_UPDATE
    %       Updates state arrays via calculated path-integral elements.
    %
    %   SYNTAX:
    %       U = path_conservative_update(U_old, I_pi, lmbd, bc, nx, nbc, inner)
    %
    %   DESCRIPTION: 
    %       Performs finite volume time update with non-conservative jump.
    %
    %   INPUTS:
    %       U_old              - Previous state vector.
    %       I_pi               - Path integral jumps.
    %       lmbd               - CFL ratio (dt/dx).
    %       bc, nx, nbc, inner - Grid/Boundary definitions.
    %
    %   OUTPUTS:
    %       U                  - Updated state vector.
    %
    %   See also: PATH_INTEGRAL_JUMPS.
    %
    %% CODE:

    U = U_old;
    for j = inner
        if bc == 2 && j == nbc+1 
            I_pi_L = I_pi(:,nx); 
        else 
            I_pi_L = I_pi(:,j-1); 
        end

        if bc == 2 && j == nx-nbc 
            I_pi_R = I_pi(:,nx); 
        else 
            I_pi_R = I_pi(:,j); 
        end

        U(:,j) = 0.5 * (U_old(:,j-1) + U_old(:,j+1)) - (lmbd/2) * (I_pi_L + I_pi_R);
    end
end