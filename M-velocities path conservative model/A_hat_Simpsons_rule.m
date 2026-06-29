function A_hat = A_hat_Simpsons_rule(U1,U2,V,N,m,R)
    %% A_HAT_SIMPSONS_RULE 
    %       Approximates path-integral matrix elements along segment 
    %       interfaces.
    %
    %   SYNTAX:
    %       A_hat = A_hat_Simpsons_rule(U1, U2, V, N, m, R)
    %
    %   INPUTS:
    %       U1, U2     - Left and right states at the interface.
    %       V, N, m, R - Numerical configuration parameters.
    %
    %   OUTPUTS:
    %       A_hat      - Integrated matrix jump approximation using 
    %                    Simpsons 1/4/1 weighting.
    %
    %   See also: A_BAR_TWO_CLASSES, PATH_INTEGRAL_JUMPS.
    %
    %% CODE:

    Umid  = 0.5*(U1+U2);
    A1    = A_bar_two_classes(U1,  V,N,m,R);
    Amid  = A_bar_two_classes(Umid,V,N,m,R);
    A2    = A_bar_two_classes(U2,  V,N,m,R);
    A_hat = (A1 + 4*Amid + A2)/6;
end