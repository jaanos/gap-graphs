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

# The polar graph O^{(+/-)}(d, q) of isotropic lines of F_q^d
# with respect to a nondegenerate quadratic form.
BindGlobal("PolarGraphO", function(arg)
    local d, e, q, B, G, H, Q, V;
    if Length(arg) < 2 then
        Error("at least two arguments expected");
        return fail;
    elif Length(arg) = 2 then
        e := 0;
        d := arg[1];
        q := arg[2];
    else
        e := arg[1];
        d := arg[2];
        q := arg[3];
    fi;
    G := GO(e, d, q);
    Q := InvariantQuadraticForm(G).matrix;
    B := Q + TransposedMat(Q);
    V := GF(q)^d;
    H := Graph(G, IsotropicSpacesQuadraticForm(Q)(Subspaces(V, 1)),
                OnSubspaces(V), function(x, y)
                    return x <> y and IsZero(Elements(x)[2]*B*Elements(y)[2]);
                end, true);
    H.duality := ProjectiveDualityFunction;
    H.primality := ProjectivePrimalityFunction;
    return H;
end);

# The polar graph Sp(d, q) of isotropic lines of F_q^d
# with respect to a nondegenerate symplectic form.
BindGlobal("PolarGraphSp", function(d, q)
    local G, H, Q, V;
    G := Sp(d, q);
    Q := InvariantBilinearForm(G).matrix;
    V := GF(q)^d;
    H := Graph(G, [Subspace(V, Elements(CanonicalBasis(V)){[1]}, "basis")],
                OnSubspaces(V), function(x, y)
                    return x <> y and IsZero(Elements(x)[2]*Q*Elements(y)[2]);
                end);
    H.duality := ProjectiveDualityFunction;
    H.primality := ProjectivePrimalityFunction;
    return H;
end);

# The polar graph U(d, r) of isotropic lines of F_{r^2}^d
# with respect to a nondegenerate Hermitean form.
BindGlobal("PolarGraphU", function(d, r)
    local c, B, F, G, H, Q, V;
    G := GU(d, r);
    Q := InvariantSesquilinearForm(G).matrix;
    V := GF(r^2)^d;
    B := Elements(CanonicalBasis(V));
    c := Conjugates(GF(r^2), GF(r), Z(r^2));
    F := x -> List(x, y -> y^r);
    H := Graph(G, [Subspace(V, [B[1] + (c[1]-c[2])*B[d]], "basis")],
            OnSubspaces(V), function(x, y)
                return x <> y and IsZero(Elements(x)[2]*Q*F(Elements(y)[2]));
            end);
    H.duality := ProjectiveDualityFunction;
    H.primality := ProjectivePrimalityFunction;
    return H;
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

# The unitary nonisotropics graph of 1-dimensional subspaces of F_(r^2)^3 with
# respect to a nondegenerate sesquilinear form.
BindGlobal("UnitaryNonisotropicsGraph", function(r)
    local G, V;
    V := GF(r^2)^3;
    G := GU(3, r);
    return Graph(G, [Subspace(V, BasisVectors(Basis(V)){[2]}, "basis")],
        OnSubspaces(V), IsOrthogonal(InvariantSesquilinearForm(G).matrix, r));
end);
