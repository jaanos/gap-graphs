RequirePackage("grape");

dpgB := function(d, q)
    local P, Q, V, i, j, e;
    e := 2*d+1;
    V := GF(q);
    Q := IdentityMat(e, V);
    if q mod 2 = 0 then
        for i in [2..e] do
            for j in [i+1..e] do
                Q[i][j] := Z(q)^0;
            od;
        od;
    fi;
    P := Filtered(Subspaces(V^e, d),
        y -> Size(Filtered(Elements(y), x -> x*Q*x <> 0*Z(q))) = 0);;
    return Graph(Group(()), P, function(x,y) return x; end,
        function(x,y)
            return Dimension(Intersection(x,y)) = d-1;
        end, true);
end;

dpgC := function(d, q)
    local I, P, Q, V, i, j, e;
    e := 2*d;
    V := GF(q);
    I := IdentityMat(d, V);
    Q := BlockMatrix([[1,2,I], [2,1,-I]], 2, 2)*Z(q)^0;
    P := Filtered(Subspaces(V^e, d),
        y -> Size(Filtered(Cartesian(Elements(y), Elements(y)),
            x -> x[1]*Q*x[2] <> 0*Z(q))) = 0);;
    return Graph(Group(()), P, function(x,y) return x; end,
        function(x,y)
            return Dimension(Intersection(x,y)) = d-1;
        end, true);
end;
