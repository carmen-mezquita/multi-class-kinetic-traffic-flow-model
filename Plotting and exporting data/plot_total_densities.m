function [pltr,pltq,ttl] = plot_total_densities(x,r,q,t,x0,xf,i)
    %% PLOT_TOTAL_DENSITIES
    %       Real-time plot updates tracking combined parameters.
    %
    %   SYNTAX:
    %       [pltr, pltq, ttl] = plot_total_densities(x, r, q, t, x0, xf, i)
    %
    %   INPUTS:
    %       x, r, q, t, x0, xf, i - Plotting configuration and data 
    %                               buffers.
    %
    %   OUTPUTS:
    %       pltr, pltq, ttl       - Handles for plot objects.
    %
    %% CODE:

    figure(i); 
    clf;
    set(gcf,'PaperPositionMode','auto'); 
    set(groot,'defaultLegendInterpreter','latex');

    subplot(2,1,1); 
    pltr = plot(x,r(:),'b'); 
    xlabel('$x$','Interpreter','latex'); 
    xlim([x0 xf]); 
    ylabel('$r$','Interpreter','latex'); 
    title('Path-conservative method');

    subplot(2,1,2); 
    pltq = plot(x,q(:),'b'); 
    xlabel('$x$','Interpreter','latex'); 
    xlim([x0 xf]); 
    ylabel('$\bar{q}^2$','Interpreter','latex'); 

    ttl = title(['time = ' num2str(t)]); 
    drawnow;
end