function func = feq(V,R,rho,r,m,nx,c)
    %% FEQ 
    %       Evaluates discrete target Maxwellian distribution vectors for 
    %       spatial relaxation steps.
    %
    %   SYNTAX:
    %       func = feq(V, R, rho, r, m, nx, c)
    %
    %   DESCRIPTION:
    %       Builds local structural equilibrium functions modeling vehicle 
    %       target velocities conditioned by spatial traffic capacity 
    %       restriction boundaries.
    %
    %   INPUTS:
    %       V, R, rho, r - Kinetic and macroscopic parameters.
    %       m, nx        - Discretization dimensions.
    %       c            - Class index.
    %
    %   OUTPUTS:
    %       func         - Target distribution vector.
    %
    %   See also: RELAXATION_TWO_CLASS.
    %
    %% CODE:

    func = zeros(m,nx);

    switch c 
        case 1 
            func(m,:) = V*rho.*max(0,(1-r/R));

        case 2
            func(m,:) = V*rho.*max(0,(1-r/R));
    end

    func(1,:) = rho - func(m,:);
end