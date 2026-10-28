function A_bar = A_bar_two_classes(U,V,N,m)
    %% A_BAR_TWO_CLASSES 
    %       Assembles system matrices for a 2-class non-conservative model.
    %
    %   SYNTAX:
    %       A_bar = A_bar_two_classes(U, V, N, m)
    %
    %   INPUTS:
    %       U          - Packed state vector.
    %       V, N, m    - Configuration and parameter inputs.
    %
    %   OUTPUTS:
    %       A_bar      - System non-conservative interaction matrix.
    %
    %   See also: AC_BLOCK, A_BAR_N_CLASSES.
    %
    %% CODE:

    f1 = U(1:m);
    f2 = U(m+1:end);
    r = sum(f2);

    Ac1 = Ac_block(f1,V,r,m);
    Ac2 = Ac_block(f2,V,r,m);

    Z = zeros(m,m);
    V_tilde = diag(repmat(V,1,N));
    A_tilde = [Z, Ac1; Z, Ac2];
    A_bar = A_tilde + V_tilde;
end