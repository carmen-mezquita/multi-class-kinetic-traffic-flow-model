function [I,total_mass,timeline,rho1,rho2,q1,q2] = discvelo_path_conservative(par)
    %% DISCVELO_PATH_CONSERVATIVE 
    %       Main execution engine for the discrete velocity traffic model.
    %
    %   SYNTAX:
    %       [I, total_mass, timeline, rho1, rho2, q1, q2] = discvelo_path_conservative(par)
    %
    %   DESCRIPTION:
    %       This core function executes a 1D, 2-class traffic simulation 
    %       based on a non-local discrete velocity kinetic framework. The 
    %       system is solved using a path-conservative finite volume scheme 
    %       to correctly evaluate non-conservative matrix jump interactions 
    %       at cell interfaces. It includes explicit steps for time-step 
    %       constraint (CFL checking), path-integral evaluation via 
    %       Simpsons rule, a conservative projection adjustment, and 
    %       non-local relaxation.
    %
    %   INPUTS:
    %       par        - Configuration struct containing domain and 
    %                    discretization definitions.
    %
    %   OUTPUTS:
    %       I          - Evolution array capturing spatial perturbation 
    %                    decay integrals.
    %       total_mass - Integrated density tracker to assess global mass 
    %                    conservation.
    %       timeline   - Array containing timestamps matching each 
    %                    simulation iteration.
    %       rho1, rho2 - Final spatial macro density tracking profiles.
    %       q1, q2     - Final momentum/flux distributions.
    %
    %   See also: PATH_INTEGRAL_JUMPS, PATH_CONSERVATIVE_UPDATE, 
    %             RELAXATION_TWO_CLASS.
    %
    %% CODE:

    % ----------- Structural Geometry Initialization -----------
    m = par.M + 1;      % Total velocity grid bins
    time = par.t_begin; % Dynamic clock register
    timeline = zeros(1,1);
    timeline(1) = time;
    it = 1;
    
    tic; % Track processing execution overhead
    
    
    % ----------- Phase-Space Distribution Allocations -----------
    [f1,f2] = deal(zeros(m,par.nx)); % Distribution matrices: [Velocity Bins x Domain Cells]
    V = linspace(0,1,m);             % Equidistant velocity configuration vector

    
    % ----------- Set Profiles and Boundary Pads -----------
    f1(:,par.inner) = initial_condition(par.x,m,V,1,par.ic); 
    f2(:,par.inner) = initial_condition(par.x,m,V,2,par.ic);
    
    f1 = boundary_conditions(par.nx,par.nbc,f1,par.bc,V);
    f2 = boundary_conditions(par.nx,par.nbc,f2,par.bc,V);

    
    % ----------- Compute Initial Macroscopic Invariants -----------
    rho1(:) = sum(f1(:,:)); 
    rho2(:) = sum(f2(:,:)); 

    global rho0
    total_mass(it) = par.dx*sum(rho2(par.inner));
    I(it) = par.dx*sum(abs(rho2(par.inner)-rho0));
    
    q1 = V*f1;
    q2 = V*f2;
    
    % Open operational file logging handles
    [pltr,pltq,ttl] = plot_total_densities(par.x_long,rho2(:),q2(:),time,par.x_low,par.x_up,1);
    file_I = fopen('I_t.txt','w');
    file_mass = fopen('total_mass.txt','w');

    
    % ----------- Main Path-Conservative Time Integration Loop -----------
    while(time<par.t_end)
        % Compute Maximum System Wave Speed for CFL Constraint
        lm = max(lambda_max(f1,V),lambda_max(f2,V));
        dt = 0.5*0.9*par.dx/lm;
        dt = min(dt,par.t_end-time);
        lmbd = dt/par.dx;
        
        f1_old = f1;
        f2_old = f2;
        U_old = [f1; f2]; % Packed matrix representation of state vectors
        
        % Step A: Evaluate path-integral matrix jumps across cellular interfaces
        I_pi = path_integral_jumps(U_old,V,par.N,m,par.bc,par.nx,par.nbc);
        
        % Step B: Apply Lax-Friedrichs type path-conservative update step
        U = path_conservative_update(U_old,I_pi,lmbd,par.bc,par.nx,par.nbc,par.inner);
        
        % Enforce basic ordering constraints on state combinations
        f1 = max(0,U(1:m,:));      % Positivity check on class 1
        f2 = max(f1,U(m+1:end,:)); % Ensure combined ordering total is >= class 1
        
        % Step C: Execute non-local kinetic relaxation step
        if par.rel==1
            [f1,f2] = relaxation_two_class(f1,f2,par.V_eps,par.R_eps,par.T_eps,m,par.nx,dt);
        end
        
        % Step D: Enforce macroscopic mass conservation properties via 
        %         projection scaling
        [f1,f2] = density_conservation([f1; f2],f1_old,f2_old,V,lmbd,par.inner,m,par.nx);
        
        % Refresh ghost boundary zones
        f1 = boundary_conditions(par.nx,par.nbc,f1,par.bc,V);
        f2 = boundary_conditions(par.nx,par.nbc,f2,par.bc,V);
        
        % Re-evaluate diagnostic variables
        rho1 = sum(f1); 
        rho2 = sum(f2);
        q1 = V*f1; 
        q2 = V*f2;
        
        time = time + dt;
        it = it + 1;
        timeline(it) = time;
        total_mass(it) = par.dx*sum(rho2(par.inner));
        I(it) = par.dx*sum(abs(rho2(par.inner)-rho0)); 
        
        data_mass = [time' total_mass(it)'];
        data_I    = [time' I(it)'];
        
        if(mod(it,4) == 0)
            fprintf(file_I,'%12.8f %12.8f\n',data_I);
            fprintf(file_mass,'%12.8f %12.8f\n',data_mass);
        end
        
        cc = 1000;
        if(mod(it,cc) == 0 || time>=par.t_end)
            pltr.YData = rho2(:);
            pltq.YData = q2(:);
            ttl.String = num2str(time);
            drawnow;
        end
    end
    CPU = toc;
    fclose(file_I);
    fclose(file_mass);
    
    plot_decay(timeline,I,total_mass,par.x_long,rho1,rho2,q1,q2,par.x_low,par.x_up,par.t_end,3);
    save_densities(par.x_long,rho1,rho2,q1,q2);
end