function A_bar = A_bar_N_classes(U,V,N,m)
    %% A_BAR_N_CLASSES 
    %       Assembles block interaction matrices generalized across N 
    %       vehicle classes.
    %
    %   SYNTAX:
    %       A_bar = A_bar_N_classes(U, V, N, m)
    %
    %   INPUTS:
    %       U          - Packed state vector.
    %       V, N, m    - Configuration and parameter inputs.
    %
    %   OUTPUTS:
    %       A_bar      - Assembled block-diagonal interaction matrix for N 
    %                    classes.
    %
    %   See also: AC_BLOCK, A_BAR_TWO_CLASSES.
    %
    %% CODE:

    fc = cell(1,N);
    for c = 1:N
        fc{c} = U((c-1)*m+1 : c*m);
    end
    r = sum(fc{N});

    A_tilde = zeros(N*m,N*m);
    A_stack = [];
    for c = 1:N
        Ac = Ac_block(fc{c},V,r,m);
        A_stack = [A_stack; Ac];
    end

    V_tilde = diag(repmat(V,1,N));
    A_tilde(:,(N-1)*m+1:end) = A_stack;
    A_bar = A_tilde + V_tilde;
end