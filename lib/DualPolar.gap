# The dual polar graph B_d(q) of the isotropic d-dimensional subspaces of
# F_q^{2d+1} with respect to a nondegenerate quadratic form.
BindGlobal("DualPolarB", function(d, q)
    local Q, V, i, j, e;
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
    return SubspaceGraph(Filtered(Subspaces(V^e, d),
        y -> Size(Filtered(Elements(y), x -> x*Q*x <> 0*Z(q))) = 0), d);
end);

# The dual polar graph C_d(q) of the isotropic d-dimensional subspaces of
# F_q^{2d} with respect to a nondegenerate symplectic form.
BindGlobal("DualPolarC", function(d, q)
    local I, Q, V, e;
    e := 2*d;
    V := GF(q);
    I := IdentityMat(d, V);
    Q := BlockMatrix([[1,2,I], [2,1,-I]], 2, 2)*Z(q)^0;
    return SubspaceGraph(Filtered(Subspaces(V^e, d),
        y -> Size(Filtered(Cartesian(Elements(y), Elements(y)),
            x -> x[1]*Q*x[2] <> 0*Z(q))) = 0), d);
end);

# The dual polar graph D_d(q) of the isotropic d-dimensional subspaces of
# F_q^{2d} with respect to a nondegenerate quadratic form of Witt index d.
BindGlobal("DualPolarD", function(d, q)
    local Q, V;
    V := GF(q);
    Q := BlockMatrix([[1, 2, IdentityMat(d, V)]], 2, 2);
    return SubspaceGraph(Filtered(Subspaces(V^(2*d), d),
        y -> Size(Filtered(Elements(y), x -> x*Q*x <> 0*Z(q))) = 0), d);
end);

# The dual polar graph ^2D_{d+1}(q) of the isotropic d-dimensional subspaces of
# F_q^{2d+2} with respect to a nondegenerate quadratic form of Witt index d.
BindGlobal("DualPolar2D", function(d, q)
    local I, Q, V, e;
    e := 2*(d+1);
    V := GF(q);
    if d = 1 then
        Q := [[0,1,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,1]]*Z(q)^0;
    else
        I := IdentityMat(d, V);
        Q := List(BlockMatrix([[1, 2, I], [3, 3, I]], 3, 3){[1..e]}, r -> r{[1..e]});
    fi;
    if q mod 2 = 0 then
        Q[e-1][e] := Z(q);
    fi;
    return SubspaceGraph(Filtered(Subspaces(V^e, d),
        y -> Size(Filtered(Elements(y), x -> x*Q*x <> 0*Z(q))) = 0), d);
end);