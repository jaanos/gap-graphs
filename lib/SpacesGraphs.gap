# The Grassmann graph J_q(n, d) of d-dimensional subspaces of F_q^n.
BindGlobal("GrassmannGraph", function(q, n, d)
    local V;
    V := GF(q)^n;
    return SubspaceGraph(GL(n, q), Elements, V, d);
end);

# The dual polar graph B_d(q) of the isotropic d-dimensional subspaces of
# F_q^{2d+1} with respect to a nondegenerate quadratic form.
BindGlobal("DualPolarGraphB", function(d, q)
    local G, V;
    V := GF(q)^(2*d+1);
    G := GO(2*d+1, q);
    return SubspaceGraph(G, IsotropicSpacesQuadraticForm(q,
        InvariantQuadraticForm(G).matrix), V, d);
end);

# The dual polar graph C_d(q) of the isotropic d-dimensional subspaces of
# F_q^{2d} with respect to a nondegenerate symplectic form.
BindGlobal("DualPolarGraphC", function(d, q)
    local G, V;
    V := GF(q)^(2*d);
    G := Sp(2*d, q);
    return SubspaceGraph(G, IsotropicSpacesBilinearForm(q,
        InvariantBilinearForm(G).matrix), V, d);
end);

# The dual polar graph D_d(q) of the isotropic d-dimensional subspaces of
# F_q^{2d} with respect to a nondegenerate quadratic form of Witt index d.
BindGlobal("DualPolarGraphD", function(d, q)
    local G, V;
    V := GF(q)^(2*d);
    G := GO(1, 2*d, q);
    return SubspaceGraph(G, IsotropicSpacesQuadraticForm(q,
        InvariantQuadraticForm(G).matrix), V, d);
end);

# The dual polar graph ^2D_{d+1}(q) of the isotropic d-dimensional subspaces of
# F_q^{2d+2} with respect to a nondegenerate quadratic form of Witt index d.
BindGlobal("DualPolarGraph2D", function(d, q)
    local G, V;
    V := GF(q)^(2*d+2);
    G := GO(-1, 2*d+2, q);
    return SubspaceGraph(G, IsotropicSpacesQuadraticForm(q,
        InvariantQuadraticForm(G).matrix), V, d);
end);
