# The Grassmann graph J_q(n, d) of d-dimensional subspaces of F_q^n.
BindGlobal("GrassmannGraph", function(q, n, d)
    return SubspaceGraph(GL(n, q), Elements, GF(q)^n, d, true);
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
    local F, e;
    e := 2*d;
    F := GF(q);
    return SubspaceGraph(Sp(e, q),
        [VectorSpace(F, Elements(CanonicalBasis(F^e)){[1..d]})],
        F^e, d, false);
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
    local B, F, c, d;
    F := GF(r^2);
    B := Elements(CanonicalBasis(F^e));
    c := Conjugates(F, GF(r), Z(r^2));
    d := Int(e/2);
    return SubspaceGraph(GU(e, r),
        [VectorSpace(F, B{[1..d]} + (c[1]-c[2])*B{e+1-[1..d]})],
        F^e, d, false);
end);
