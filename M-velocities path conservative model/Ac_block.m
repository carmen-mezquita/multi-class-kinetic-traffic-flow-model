function Ac = Ac_block(fc,V,r,m)
    %% AC_BLOCK 
    %       Evaluates directional non-conservative interaction derivatives 
    %       inside cells.
    %
    %   SYNTAX:
    %       Ac = Ac_block(fc, V, r, m)
    %
    %   DESCRIPTION:
    %       Assembles lower-triangular velocity jump blocks handling 
    %       non-local velocity adjustments
    %
    %   INPUTS:
    %       fc    - Class distribution.
    %       V     - Velocity vector.
    %       r     - Total density.
    %       m     - Velocity bins and capacity limit.
    %
    %   OUTPUTS:
    %       Ac    - Lower-triangular velocity jump block matrix.
    %
    %   See also: A_BAR_TWO_CLASSES, A_BAR_N_CLASSES.
    %
    %% CODE:

    Ac = zeros(m,m);
    if r<1
        denom = max(1-r,1e-12);
        for i = 1:m
            cte = 0;
            for l = (i+1):m
                cte = cte + (V(l)-V(i))*fc(l);
            end
            Ac(i,i) = -cte/denom;
            if i > 1
                Ac(i,1:i-1) = (V(i)-V(1:i-1))*fc(i)/denom;
            end
        end
    end
end