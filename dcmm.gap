RequirePackage("grape");

dCMM := function(t, h, e)
    local q, s, F, K;
    q := 2^(2*t-1);
    s := 2^e;
    F := GF(q);
    if h = 1 then
        return Graph(Group(()), Cartesian(F, GF(2), F),
            function(x,y) return x; end,
            function(x,y)
                return x <> y
                    and x[3]+y[3] = x[1]^s * y[1] + x[1] * y[1]^s + (x[2]+y[2])*(x[1]^(s+1) + y[1]^(s+1));
            end, true);
    else
        K := AdditiveGroup(BasisVectors(Basis(F)){[1..h-1]});
        return Graph(Group(()), Cartesian(F, GF(2), Unique(List(F, x -> x+K))),
            function(x,y) return x; end,
            function(x,y)
                return x <> y
                    and x[3]+Elements(y[3])[1] + x[1]^s * y[1] + x[1] * y[1]^s + (x[2]+y[2])*(x[1]^(s+1) + y[1]^(s+1)) = K;
            end, true);
    fi;
end;