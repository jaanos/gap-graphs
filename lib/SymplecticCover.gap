RequirePackage("grape");

AdditiveSympCover := function(q, n, h)
    local B, F, K, V, i, j;
    F := GF(q);
    V := F^n;
    B := InvariantBilinearForm(Sp(n, q)).matrix;
    if h = 0 then
        return Graph(Group(()), Cartesian(F, V),
            function(x, y) return x; end,
            function(x, y)
                return x <> y and x[2]*B*y[2]+y[1]-x[1] = 0*Z(q);
            end, true);
    else
        K := AdditiveGroup(BasisVectors(Basis(F)){[1..h]});
        return Graph(Group(()), Cartesian(Unique(List(F, x -> x+K)), V),
            function(x, y) return x; end,
            function(x, y)
                return x <> y and x[2]*B*y[2]+y[1]-Elements(x[1])[1] = K;
            end, true);
    fi;
end;

MultiplicativeSympCover := function(q, m)
    local B, F, K;
    F := GF(q);
    K := Elements(Group(Z(q)^((q-1)/m)));
    B := [[0,1],[-1,0]]*Z(q)^0;
    return Graph(Group(()), Unique(List(Filtered(F^2, x -> x <> Zero(F^2)),
            v -> Set(List(K, g -> List(v, x -> g*x))))),
        function(x, y) return x; end,
        function(x, y)
            return x <> y and x[1]*B*y[1] in K;
        end, true);
end;