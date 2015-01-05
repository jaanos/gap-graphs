# The Grassmann graph J_q(n, d) of d-dimensional subspaces of F_q^n.
BindGlobal("GrassmannGraph", function(q, n, d)
    return SubspaceGraph(GL(n, q), Elements, GF(q)^n, d, true);
end);

# The twisted Grassmann graph TG_d(q) of (d+1)-dimensional subspaces of
# F_q^{2d+1} which are not subspaces of a hyperplane H, and (d-1)-dimensional
# subspaces of H.
BindGlobal("TwistedGrassmannGraph", function(d, q)
    local F, H, S, V, n;
    n := 2*d+1;
    F := GF(q)^n;
    H := Subspace(F, CanonicalBasis(F){[1..n-1]}, "basis");
    V := Union(Filtered(Subspaces(F, d+1), x -> not IsSubset(H, x)),
            Subspaces(H, d-1));
    S := OnSubspaces(F);
    return Graph(Stabilizer(GL(n, q), H, S), V, S, function(x, y)
        local d1, d2;
        d1 := Dimension(x);
        d2 := Dimension(y);
        if d1 < d2 then
            return IsSubset(y, x);
        elif d1 > d2 then
            return IsSubset(x, y);
        else
            return Dimension(Intersection(x, y)) = d1-1;
        fi;
    end, true);
end);

# The dual polar graph B_d(q) of isotropic d-dimensional subspaces of
# F_q^{2d+1} with respect to a nondegenerate quadratic form.
BindGlobal("DualPolarGraphB", function(d, q)
    local G, e;
    e := 2*d+1;
    G := GO(e, q);
    return SubspaceGraph(G,
        IsotropicSpacesQuadraticForm(InvariantQuadraticForm(G).matrix),
        GF(q)^e, d);
end);

# The dual polar graph C_d(q) of isotropic d-dimensional subspaces of
# F_q^{2d} with respect to a nondegenerate symplectic form.
BindGlobal("DualPolarGraphC", function(d, q)
    local F, e;
    e := 2*d;
    F := GF(q)^e;
    return SubspaceGraph(Sp(e, q),
        [Subspace(F, Elements(CanonicalBasis(F)){[1..d]}, "basis")],
        F, d, false);
end);

# The dual polar graph D_d(q) of isotropic d-dimensional subspaces of
# F_q^{2d} with respect to a nondegenerate quadratic form of Witt index d.
BindGlobal("DualPolarGraphD", function(d, q)
    local G, e;
    e := 2*d;
    G := GO(1, e, q);
    return SubspaceGraph(G,
        IsotropicSpacesQuadraticForm(InvariantQuadraticForm(G).matrix),
        GF(q)^e, d);
end);

# The dual polar graph ^2D_{d+1}(q) of isotropic d-dimensional subspaces of
# F_q^{2d+2} with respect to a nondegenerate quadratic form of Witt index d.
BindGlobal("DualPolarGraph2D", function(d, q)
    local G, e;
    e := 2*d+2;
    G := GO(-1, e, q);
    return SubspaceGraph(G,
        IsotropicSpacesQuadraticForm(InvariantQuadraticForm(G).matrix),
        GF(q)^e, d);
end);

# The dual polar graph ^2A_{e-1}(r) of isotropic [e/2]-dimensional
# subspaces of F_{r^2}^e with respect to a nondegenerate Hermitean form.
BindGlobal("DualPolarGraph2A", function(e, r)
    local B, F, c, d;
    F := GF(r^2)^e;
    B := Elements(CanonicalBasis(F));
    c := Conjugates(GF(r^2), GF(r), Z(r^2));
    d := Int(e/2);
    return SubspaceGraph(GU(e, r),
        [Subspace(F, B{[1..d]} + (c[1]-c[2])*B{e+1-[1..d]}, "basis")],
        F, d, false);
end);

# The Doro graph of nonisotropic 1-dimensional subspaces of F_q^4 with respect
# to a nondegenerate quadratic form. It is distance-regular for q = 4, 5.
BindGlobal("DoroGraph", function(q)
    local G, V;
    V := GF(q)^4;
    G := GO(-1, 4, q);
    return Graph(G, [Subspace(V, BasisVectors(Basis(V)){[4]}, "basis")],
        OnSubspaces(V), IsHyperbolic(InvariantQuadraticForm(G).matrix));
end);

