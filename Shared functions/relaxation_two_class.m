function [f1,f2] = relaxation_two_class(f1,f2,V_eps,R_eps,T_eps,m,nx,dt)
    %% RELAXATION_TWO_CLASS 
    %       Executes non-local interaction updates via relaxation time 
    %       scales.
    %
    %   SYNTAX:
    %       [f1, f2] = relaxation_two_class(f1, f2, V_eps, R_eps, T_eps, m, nx, dt)
    %
    %   DESCRIPTION:
    %       Maps modified system distributions to physical configurations, 
    %       evaluates local equilibrium targets via 'feq', solves source 
    %       terms over timestamp dt, and re-maps back.
    %
    %   INPUTS:
    %       f1, f2              - Current class distributions.
    %       V_eps, R_eps, T_eps - Relaxation parameters.
    %       m, nx, dt           - Discretization and time parameters.
    %
    %   OUTPUTS:
    %       f1, f2              - Updated (relaxed) distributions.
    %
    %   See also: FEQ.
    %
    %% CODE:

    rho1 = sum(f1); 
    rho2 = sum(f2); 

    f1aux = f1;           
    f2aux = max(0,f2-f1); 
    
    rho1aux = rho1;        
    rho2aux = rho2 - rho1; 
    
    feq1 = feq(V_eps(1),R_eps(1),rho1aux,rho2,m,nx,1);
    feq2 = feq(V_eps(2),R_eps(2),rho2aux,rho2,m,nx,2);
    
    f1aux = (T_eps(1)*f1aux + dt*feq1)/(T_eps(1) + dt);
    f2aux = (T_eps(2)*f2aux + dt*feq2)/(T_eps(2) + dt);
    
    f1 = f1aux;
    f2 = f1aux + f2aux;
end