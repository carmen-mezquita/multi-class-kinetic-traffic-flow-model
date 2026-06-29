function plot_decay(timeline,I,total_mass,x,rho1,rho2,q1,q2,x0,xf,tf,i)
    %% PLOT_DECAY 
    %       Exports finalized diagnostics tracking perturbation decay 
    %       invariants.
    %
    %   SYNTAX:
    %       plot_decay(timeline, I, total_mass, x, rho1, rho2, q1, q2, x0, xf, tf, i)
    %
    %   INPUTS:
    %       timeline, I, total_mass - Time-evolution data.
    %       x, rho1, rho2, q1, q2   - Final state arrays.
    %       x0, xf, tf, i           - Plotting and domain constraints.
    %
    %% CODE:

    global rho0

    figI = figure(i); 
    plot(timeline,I,'LineWidth',2); 
    ytxt = sprintf('$I(t)=\\int | \\rho(x,t) - %g | \\, dx$', rho0);
    xlabel('$t$','Interpreter','latex'); 
    ylabel(ytxt,'Interpreter','latex'); 
    title('Decay of density perturbation (Test-case 2)'); 
    exportgraphics(figI,'I_t.png','Resolution',300);
    
    fig_mass = figure(i+1); 
    plot(timeline,total_mass,'LineWidth',2); 
    xlabel('$t$','Interpreter','latex'); 
    ylabel('$r$','Interpreter','latex'); 
    title('Total mass (Test-case 2)'); 
    exportgraphics(fig_mass,'total_mass.png','Resolution',300);
    
    fig_class_densities = figure(i+2); 
    clf; 
    set(gcf,'PaperPositionMode','auto'); 
    set(groot,'defaultLegendInterpreter','latex');

    subplot(2,1,1); 
    plot(x,rho1(:),'r',x,rho2(:)-rho1(:),'b'); 
    xlabel('$x$','Interpreter','latex'); 
    xlim([x0 xf]); 
    ylabel('$\rho^c$','Interpreter','latex'); 
    leg1 = legend('$\rho^1$','$\rho^2$'); 
    set(leg1,'Interpreter','latex'); 
    title('Path-conservative method');

    subplot(2,1,2); 
    plot(x,q1(:),'r',x,q2(:)-q1(:),'b'); 
    xlabel('$x$','Interpreter','latex'); 
    xlim([x0 xf]); 
    ylabel('$q^c$','Interpreter','latex'); 
    leg2 = legend('$q^1$','$q^2$'); 
    set(leg2,'Interpreter','latex'); 
    title(['time = ' num2str(tf)]); 
    drawnow;
    
    exportgraphics(fig_class_densities,'class_densities.png','Resolution',300);
end