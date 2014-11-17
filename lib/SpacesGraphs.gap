# The Grassmann graph J_q(n, d) of d-dimensional subspaces of F_q^n.
BindGlobal("GrassmannGraph", function(q, n, d)
    return SubspaceGraph(GL(n, q), Elements, GF(q)^n, d);
end);

# The dual polar graph B_d(q) of the isotropic d-dimensional subspaces of
# F_q^{2d+1} with respect to a nondegenerate quadratic form.
BindGlobal("DualPolarGraphB", function(d, q)
    local G, e;
    e := 2*d+1;
    G := GO(e, q);
    return SubspaceGraph(G,
        IsotropicSpacesQuadraticForm(InvariantQuadraticForm(G).matrix),
        GF(q)^e, d);
end);

# The dual polar graph C_d(q) of the isotropic d-dimensional subspaces of
# F_q^{2d} with respect to a nondegenerate symplectic form.
BindGlobal("DualPolarGraphC", function(d, q)
    local G, e;
    e := 2*d
    G := Sp(e, q);
    return SubspaceGraph(G,
        IsotropicSpacesBilinearForm(InvariantBilinearForm(G).matrix),
        GF(q)^e, d);
end);

# The dual polar graph D_d(q) of the isotropic d-dimensional subspaces of
# F_q^{2d} with respect to a nondegenerate quadratic form of Witt index d.
BindGlobal("DualPolarGraphD", function(d, q)
    local G, e;
    e := 2*d;
    G := GO(1, e, q);
    return SubspaceGraph(G,
        IsotropicSpacesQuadraticForm(InvariantQuadraticForm(G).matrix),
        GF(q)^e, d);
end);

# The dual polar graph ^2D_{d+1}(q) of the isotropic d-dimensional subspaces of
# F_q^{2d+2} with respect to a nondegenerate quadratic form of Witt index d.
BindGlobal("DualPolarGraph2D", function(d, q)
    local G, e;
    e := 2*d+2;
    G := GO(-1, e, q);
    return SubspaceGraph(G,
        IsotropicSpacesQuadraticForm(InvariantQuadraticForm(G).matrix),
        GF(q)^e, d);
end);

# The dual polar graph ^2A_{e-1}(r) of the isotropic [e/2]-dimensional
# subspaces of F_{r^2}^e with respect to a nondegenerate Hermitean form.
BindGlobal("DualPolarGraph2A", function(e, r)
    local F, G;
    F := GF(r^2);
    G := GU(e, r);
    return SubspaceGraph(G,
        IsotropicSpacesSesquilinearForm(InvariantSesquilinearForm(G).matrix, r),
        F^e, Int(e/2));
end);
