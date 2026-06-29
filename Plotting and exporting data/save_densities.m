function save_densities(x,rho1,rho2,q1,q2)
    %% SAVE_DENSITIES 
    %       Formats internal numerical profiles to plain text arrays.
    %
    %   SYNTAX:
    %       save_densities(x, rho1, rho2, q1, q2)
    %
    %   INPUTS:
    %       x                  - Spatial domain vector.
    %       rho1, rho2, q1, q2 - State variable arrays to be saved.
    %
    %% CODE:

    rho_classes = [x' rho1(:) (rho2(:)-rho1(:))];
    save('rho_classes.txt','rho_classes','-ascii');

    q_classes = [x' q1(:) (q2(:)-q1(:))];
    save('q_classes.txt','q_classes','-ascii');

    r_total = [x' rho2(:)];
    save('r_total.txt','r_total','-ascii');

    q_total = [x' q2(:)];
    save('q_total.txt','q_total','-ascii');
end