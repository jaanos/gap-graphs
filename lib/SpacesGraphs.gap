# The Grassmann graph J_q(n, d) of d-dimensional subspaces of F_q^n.
BindGlobal("GrassmannGraph", function(q, n, d)
    local V;
    V := GF(q)^n;
    return Graph(GL(n, q), Elements(Subspaces(V, d)), OnSubspaces(V),
        function(x,y)
            return Dimension(Intersection(x,y)) = d-1;
        end, true);
end);

# The dual polar graph B_d(q) of the isotropic d-dimensional subspaces of
# F_q^{2d+1} with respect to a nondegenerate quadratic form.
BindGlobal("DualPolarGraphB", function(d, q)
    local Q, V, e;
    e := 2*d+1;
    V := GF(q);
    if q mod 2 = 0 then
        Q := Concatenation([Concatenation([Z(q)^0],
            ListWithIdenticalEntries(e-1, 0*Z(q)))], List([1..e-1],
                i -> Concatenation(ListWithIdenticalEntries(i, 0*Z(q)),
                                ListWithIdenticalEntries(e-i, Z(q)^0))));
    else
        Q := IdentityMat(e, V);
    fi;
    return SubspaceGraph(Filtered(Subspaces(V^e, d),
        y -> Size(Filtered(Elements(y), x -> x*Q*x <> 0*Z(q))) = 0), d);
end);

# The dual polar graph C_d(q) of the isotropic d-dimensional subspaces of
# F_q^{2d} with respect to a nondegenerate symplectic form.
BindGlobal("DualPolarGraphC", function(d, q)
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
BindGlobal("DualPolarGraphD", function(d, q)
    local Q, V;
    V := GF(q);
    Q := BlockMatrix([[1, 2, IdentityMat(d, V)]], 2, 2);
    return SubspaceGraph(Filtered(Subspaces(V^(2*d), d),
        y -> Size(Filtered(Elements(y), x -> x*Q*x <> 0*Z(q))) = 0), d);
end);

# The dual polar graph ^2D_{d+1}(q) of the isotropic d-dimensional subspaces of
# F_q^{2d+2} with respect to a nondegenerate quadratic form of Witt index d.
BindGlobal("DualPolarGraph2D", function(d, q)
    local I, J, Q, V, e;
    e := 2*(d+1);
    V := GF(q);
    I := IdentityMat(2, V);
    if q mod 2 = 0 then
        I[2][1] := Z(q);
    fi;
    J := [[0, 1], [0, 0]]*Z(q)^0;
    Q := BlockMatrix(Concatenation(List([1..d], i -> [i, i, J]),
        [[d+1, d+1, I]]), d+1, d+1);
    return SubspaceGraph(Filtered(Subspaces(V^e, d),
        y -> Size(Filtered(Elements(y), x -> x*Q*x <> 0*Z(q))) = 0), d);
end);
