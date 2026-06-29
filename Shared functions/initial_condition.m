function f_0 = initial_condition(x,m,V,c,ic)
    %% INITIAL_CONDITION 
    %       Sets spatial-density phase profiles across specific benchmark 
    %       cases.
    %
    %   SYNTAX:
    %       f_0 = initial_condition(x, m, V, c, ic)
    %
    %   DESCRIPTION:
    %       Generates initial velocity distribution profiles f_i(x,0) 
    %       depending on different initial conditions.
    %
    %   INPUTS:
    %       x     - Spatial coordinate array.
    %       m     - Velocity grid bins.
    %       V     - Velocity vector.
    %       c     - Vehicle class index.
    %       ic    - Initial condition case index.
    %
    %   OUTPUTS:
    %       f_0   - Initial distribution profile.
    %
    %   See also: DISCVELO_PATH_CONSERVATIVE.
    %
    %% CODE:

    global rho0
    rho0 = 0;
    nx_inner = length(x);
    f_0 = zeros(m,nx_inner);
    
    switch ic
        case 1 % Riemann problem, no relaxation (test case 1)
            if c==1
                rhoL = 0.4;
                rhoR = 0.2;
            elseif c==2
                rhoL = 0.6;
                rhoR = 0.8; 
            end

            FL = rhoL/(m)*ones(m,nx_inner);
            FR = rhoR/(m)*ones(m,nx_inner);
            f_0(:,:) = FL(:,:);
            f_0(:,floor(nx_inner/2):end) = FR(:,floor(nx_inner/2):end);

        case 2 % Riemann problem, relaxation (test case 2, overtaking)
            rho = zeros(1,length(x));
            if c==1
                rho = 0.5*(x>=1 & x<=10);

            elseif c==2
                rho = 0.5*((x>=1 & x<=10) + (x>=11 & x<=20));
            end

            rho0 = 0;
            for i = 1:m
                f_0(i,:) = rho/m;
            end

        case 3 % Periodic conditions, relaxation (test case 3)
            if c==1 
                rho = 0.3 + 0.1*sin(6*pi*x); 

            elseif c==2 
                rho0 = 0.7;
                rho = 0.7 + 0.25*sin(6*pi*x); 
            end

            for i = 1:m
                f_0(i,:) = rho/m;
            end  
    end
end