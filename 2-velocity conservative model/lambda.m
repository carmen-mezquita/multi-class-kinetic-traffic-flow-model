function func = lambda(u,V,m,nx)
    %% LAMBDA  
    %       Computes the characteristic velocities (eigenvalues) for the
    %       conservative formulation.
    %
    %   SYNTAX:
    %       func = lambda(u, V, m, nx)
    %
    %   DESCRIPTION: 
    %       Computes spectral speeds for the cumulative state system.
    %
    %   INPUTS:
    %       u      - Cumulative state variables.
    %       V      - Velocity vector.
    %       m, nx  - Dimensions.
    %
    %   OUTPUTS:
    %       func   - Wave speed matrix.
    %
    % See also: FSYS

    %% CODE:
    
    summation = zeros(m,nx-1);
    spd = zeros(m,nx-1);

    for i = 1:m-1
        for j = i+1:m-1
           summation(i,:) = summation(i,:) ...
                            + (V(j) - V(i)).*(1./u(j,:) - 1./u(j+1,:));
        end

        summation(i,:) = summation(i,:) + (1 - V(i)).*(1./u(m,:) - 1); 
        spd(i,:) = V(i) - summation(i,:);
    end
    
    spd(m,:) = ones(1,nx-1); % lambda_m = 1;
    func = spd;
end