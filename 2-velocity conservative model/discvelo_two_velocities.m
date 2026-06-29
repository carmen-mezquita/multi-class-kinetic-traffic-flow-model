function [I,total_mass,timeline,rho1,rho2,q1,q2] = discvelo_two_velocities(par)
    %% DISCVELO_TWO_VELOCITIES 
    %       Main execution engine for the conservative 2-velocity model.
    %
    %   SYNTAX:
    %       [I, total_mass, timeline, x, rho1, rho2, q1, q2] = discvelo_two_velocities(par)
    %
    %   DESCRIPTION:
    %       Simulates 1D, 2-class traffic flow using a strictly 
    %       conservative discrete velocity approach. The system solves for 
    %       the spatial evolution of macroscopic quantities and 
    %       distribution variables using a FORCE finite volume scheme (a 
    %       convex combination of Lax-Friedrichs and Lax-Wendroff numerical 
    %       fluxes).
    %
    %   INPUTS:
    %       par        - Unified configuration struct containing domain, 
    %                    physics, and grid definitions.
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
    %   See also: FLUX_FORCE, RECONSTRUCT_F1, RECONSTRUCT_F2, 
    %             RELAXATION_TWO_CLASS.
    %
    %% CODE:

    m = par.M + 1; 
    time = par.t_begin; 
    timeline = zeros(1,1);
    timeline(1) = time;
    it = 1;

    tic;

    % --- Phase-Space Configuration ---
    V = linspace(0,1,m);
    [u1,u2,f1,f2] = deal(zeros(m,par.nx));

    % --- Set Initial Conditions ---
    f1(:,par.inner) = initial_condition(par.x,m,V,1,par.ic);
    f2(:,par.inner) = initial_condition(par.x,m,V,2,par.ic);
    
    % Note: Passes 'V' to match the shared boundary_conditions signature
    f1 = boundary_conditions(par.nx,par.nbc,f1,par.bc,V);
    f2 = boundary_conditions(par.nx,par.nbc,f2,par.bc,V);
    
    % --- Establish Base Conservative Invariants ---
    rho1(:) = f1(1,:) + f1(2,:);
    rho2(:) = f2(1,:) + f2(2,:);

    global rho0
    total_mass(it) = par.dx*sum(rho2(par.inner));
    I(it) = par.dx*sum(abs(rho2(par.inner)-rho0));
    
    q1 = V*f1;
    q2 = V*f2; 

    % Evaluate macroscopic variable structures for both classes
    u1(1,:) = rho1(:);
    u2(1,:) = 1 - rho2(:);
    
    u1(2,:) = u2(1,:)./f1(2,:);
    u2(2,:) = u2(1,:)./(1-f2(1,:));
    
    [pltr,pltq,ttl] = plot_total_densities(par.x_long,rho2(:),q2(:),time,par.x_low,par.x_up,1);

    file_I = fopen('I_t.txt','w');
    file_mass = fopen('total_mass.txt','w');

    % --- Main FORCE Finite Volume Integration Loop ---
    while(time<par.t_end)
        % CFL
        lm1 = lambda_max(f1,V);
        lm2 = lambda_max(f2,V);
        lm = max(lm1,lm2);
        dt = 0.5*0.9*par.dx/lm;
        dt = min(dt,par.t_end-time);
        lmbd = dt/par.dx;
        
        % Evaluate interface numerical fluxes
        [flux1,flux2] = flux_FORCE(u1(:,par.inter),u1(:,par.inter+1),...
                                   u2(:,par.inter),u2(:,par.inter+1),lmbd,V,m,par.nx);

        % Apply flux gradients to update the internal cumulative parameters
        u1(:,par.inner) = u1(:,par.inner) - lmbd*(flux1(:,2:end) - flux1(:,1:end-1));
        u2(:,par.inner) = u2(:,par.inner) - lmbd*(flux2(:,2:end) - flux2(:,1:end-1));

        f1_old = f1;
        f2_old = f2;
        alpha = 1;

        % Reconstruct density arrays backwards from the cumulative variables
        [f1_target,~] = reconstruct_f1(u1,u2,m,par.nx,par.nbc,par.bc);
        [f2_target,~] = reconstruct_f2(u2,m,par.nx,par.nbc,par.bc);

        % Smooth the distribution updates
        f1 = (1 - alpha) * f1_old + alpha * f1_target;
        f2 = (1 - alpha) * f2_old + alpha * f2_target;

        rho1 = sum(f1(:,:));
        rho2 = sum(f2(:,:));

        % Enforce non-local kinetic relaxation effects
        if par.rel==1
            [f1,f2] = relaxation_two_class(f1,f2,par.V_eps,par.R_eps,par.T_eps,m,par.nx,dt);

            rho1(:) = sum(f1(:,:));
            rho2(:) = sum(f2(:,:));

            % Rectify cumulative representations to match new relaxed 
            % distributions
            u1(1,:) = rho1; 
            u2(1,:) = 1 - rho2;
            u1(2,:) = u2(1,:)./f1(2,:);
            u2(2,:) = u2(1,:)./(1-f2(1,:));
        end

        u1 = boundary_conditions(par.nx,par.nbc,u1,par.bc,V);
        u2 = boundary_conditions(par.nx,par.nbc,u2,par.bc,V);

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