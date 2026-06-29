function lammax = lambda_max(f,V)
    %% LAMBDA_MAX 
    %       Extracts largest transport/characteristic wave velocity 
    %       component.
    %
    %   SYNTAX:
    %       lammax = lambda_max(f, V)
    %
    %   DESCRIPTION:
    %       Computes the upper bound for characteristic wave speeds across
    %       the discrete grid cells to determine structural dynamic CFL 
    %       constraints.
    %
    %   INPUTS:
    %       f      - Distribution matrix.
    %       V      - Velocity vector.
    %
    %   OUTPUTS:
    %       lammax - Calculated maximum wave speed for CFL constraints.
    %
    %   See also: DISCVELO_PATH_CONSERVATIVE.
    %
    %% CODE:

    rho = sum(f);
    q = V*f;
    spd = -q./(1-rho);
    lammax = max(abs(spd));
    lammax = max(1,lammax); % Guarantee minimum lower clipping bounds
end