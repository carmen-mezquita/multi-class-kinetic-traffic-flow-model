function I_pi = path_integral_jumps(U,V,N,m,bc,nx,nbc)
    %% PATH_INTEGRAL_JUMPS 
    %       Computes path conservative matrix jumps at interfaces.
    %
    %   SYNTAX:
    %       I_pi = path_integral_jumps(U, V, N, m, bc, nx, nbc)
    %
    %   DESCRIPTION: 
    %       Integrates the non-conservative path across cell interfaces.
    %
    %   INPUTS:
    %       U                           - State vector.
    %       V, N, m, bc, nx, nbc        - Numerical parameters.
    %
    %   OUTPUTS:
    %       I_pi                        - Matrix jump approximations.
    %
    %   See also: A_HAT_SIMPSONS_RULE, PATH_CONSERVATIVE_UPDATE.
    %
    %% CODE:

    I_pi = zeros(N*m,nx); 
    for j = 1:(nx-1)
        A_hat = A_hat_Simpsons_rule(U(:,j),U(:,j+1),V,N,m);
        I_pi(:,j) = A_hat*(U(:,j+1)-U(:,j));
    end

    if bc == 2  % Periodic boundary conditions
        A_hat = A_hat_Simpsons_rule(U(:,nx-nbc),U(:,nbc+1),V,N,m);
        I_pi(:,nx) = A_hat * (U(:,nbc+1) - U(:,nx-nbc));
    end
end