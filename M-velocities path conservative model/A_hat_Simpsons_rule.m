function A_hat = A_hat_Simpsons_rule(U1,U2,V,N,m)
    %% A_HAT_SIMPSONS_RULE 
    %       Approximates path-integral matrix elements along segment 
    %       interfaces.
    %
    %   SYNTAX:
    %       A_hat = A_hat_Simpsons_rule(U1, U2, V, N, m)
    %
    %   INPUTS:
    %       U1, U2     - Left and right states at the interface.
    %       V, N, m    - Numerical configuration parameters.
    %
    %   OUTPUTS:
    %       A_hat      - Integrated matrix jump approximation using 
    %                    Simpsons 1/4/1 weighting.
    %
    %   See also: A_BAR_TWO_CLASSES, PATH_INTEGRAL_JUMPS.
    %
    %% CODE:

    d1 = 1e-8;
    d2 = 1e-12;
    U1e = U1;  
    U2e = U2;
    
    if (sum(U1(m+1:end)) >= 1 - d2)
        U1e = U1 + d1*(U2-U1); 
    end
    
    if (sum(U2(m+1:end)) >= 1 - d2)
        U2e = U2 - d1*(U2-U1); 
    end

    Umid  = 0.5*(U1+U2);
    A1    = A_bar_two_classes(U1e, V,N,m);
    Amid  = A_bar_two_classes(Umid,V,N,m);
    A2    = A_bar_two_classes(U2e, V,N,m);
    A_hat = (A1 + 4*Amid + A2)/6;
end