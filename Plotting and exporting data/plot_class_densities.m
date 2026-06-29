function [pltr,pltq,ttl] = plot_class_densities(x,rho1,rho2,q1,q2,t,x0,xf,i)
    %% PLOT_CLASS_DENSITIES 
    %       Plot updates tracking individual vehicle classes independently.
    %
    %   SYNTAX:
    %       [pltr, pltq, ttl] = plot_class_densities(x, rho1, rho2, q1, q2, t, x0, xf, i)
    %
    %   INPUTS:
    %       x, rho1, rho2, q1, q2, t, x0, xf, i - Data and formatting 
    %                                             parameters.
    %
    %   OUTPUTS:
    %       pltr, pltq, ttl                     - Handles for plot objects.
    %
    %% CODE:

    figure(i); 
    clf;
    set(gcf,'PaperPositionMode','auto'); 
    set(groot,'defaultLegendInterpreter','latex');

    subplot(2,1,1); 
    pltr = plot(x,rho1(:),'r',x,rho2(:)-rho1(:),'b'); 
    xlabel('$x$','Interpreter','latex'); 
    xlim([x0 xf]); 
    ylabel('$\rho^c$','Interpreter','latex'); 
    leg1 = legend('$\rho^1$','$\rho^2$'); 
    set(leg1,'Interpreter','latex'); 
    
    title('Path-conservative method');

    subplot(2,1,2); 
    pltq = plot(x,q1(:),'r',x,q2(:)-q1(:),'b'); 
    xlabel('$x$','Interpreter','latex'); 
    xlim([x0 xf]); 
    ylabel('$q^c$','Interpreter','latex'); 
    leg2 = legend('$q^1$','$q^2$'); 
    set(leg2,'Interpreter','latex'); 
    ttl = title(['time = ' num2str(t)]); 
    drawnow;
end