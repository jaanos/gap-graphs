# The Grassmann graph J_q(n, d) of d-dimensional subspaces of F_q^n.
InstallMethod(GrassmannGraphCons,
    "as a spaces graph with full automorphism group", true,
    [IsSpacesGraph and FullAutomorphismGroup, IsInt, IsInt, IsInt], 0,
    function(filter, q, n, d)
        local G, V, m, dp, tr;
        V := GF(q)^n;
        if d in [1, n-1] then
            m := (q^n-1)/(q-1);
            G := CompleteGraph(SymmetricGroup(m), m);
            AssignVertexNames(G, Elements(Subspaces(V, d)));
        else
            if n = 2*d then
                tr := (1, 2);
            else
                tr := ();
            fi;
            dp := DirectProduct(GL(n, q), Group(tr),
                                FieldExponentiationPermutationGroup(q));
            G := SubspaceGraph(dp, Elements, V, d, OnGrassmann(q, V, dp),
                                true);
        fi;
        G.duality := GrassmannDualityFunction;
        G.primality := GrassmannDualityFunction;
        return G;
    end);

InstallMethod(GrassmannGraphCons, "as a spaces graph", true,
    [IsSpacesGraph, IsInt, IsInt, IsInt], 0, function(filter, q, n, d)
        local G;
        G := SubspaceGraph(GL(n, q), Elements, GF(q)^n, d, true);
        G.duality := GrassmannDualityFunction;
        G.primality := GrassmannDualityFunction;
        return G;
    end);

InstallMethod(GrassmannGraphCons, "with full automorphism group", true,
    [FullAutomorphismGroup, IsInt, IsInt, IsInt], 0, function(filter, q, n, d)
        return GrassmannGraphCons(IsSpacesGraph and FullAutomorphismGroup,
                                    q, n, d);
    end);

InstallMethod(GrassmannGraphCons, "default", true,
    [IsObject, IsInt, IsInt, IsInt], 0, function(filter, q, n, d)
        return GrassmannGraphCons(IsSpacesGraph, q, n, d);
    end);

BindGlobal("GrassmannGraph", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j+2 then
        return GrassmannGraphCons(filt, arg[j], arg[j+1], arg[j+2]);
    else
        Error("usage: GrassmannGraph( [<filter>, ]<int>, <int>, <int> )");
    fi;
end);

# The doubled Grassmann graph 2J_q(2d+1, d) of d- and (d+1)-dimensional
# subspaces of F_q^{2d+1}.
InstallMethod(DoubledGrassmannGraphCons,
    "as a spaces graph with full automorphism group", true,
    [IsSpacesGraph and FullAutomorphismGroup, IsInt, IsInt], 0,
    function(filter, q, d)
        local n, G, V, dp;
        n := 2*d+1;
        V := GF(q)^n;
        dp := DirectProduct(GL(n, q), Group((1, 2)),
                            FieldExponentiationPermutationGroup(q));
        G := Graph(dp, Union(Subspaces(V, d), Subspaces(V, d+1)),
                    OnGrassmann(q, V, dp), SymmetrizedInclusion, true);
        G.halfDuality := GrassmannDualityFunction;
        G.halfPrimality := GrassmannDualityFunction;
        return G;
    end);

InstallMethod(DoubledGrassmannGraphCons, "as a spaces graph", true,
    [IsSpacesGraph, IsInt, IsInt], 0, function(filter, q, d)
        return DoubledGrassmannGraphCons(IsSpacesGraph
                                            and FullAutomorphismGroup, q, d);
    end);

InstallMethod(DoubledGrassmannGraphCons, "with full automorphism group", true,
    [FullAutomorphismGroup, IsInt, IsInt], 0, function(filter, q, d)
        return DoubledGrassmannGraphCons(IsSpacesGraph
                                            and FullAutomorphismGroup, q, d);
    end);

InstallMethod(DoubledGrassmannGraphCons, "default", true,
    [IsObject, IsInt, IsInt], 0, function(filter, q, d)
        return DoubledGrassmannGraphCons(IsSpacesGraph, q, d);
    end);

BindGlobal("DoubledGrassmannGraph", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j+1 then
        return DoubledGrassmannGraphCons(filt, arg[j], arg[j+1]);
    else
        Error("usage: DoubledGrassmannGraph( [<filter>, ]<int>, <int>, <int> )");
    fi;
end);

# The twisted Grassmann graph TG_d(q) of (d+1)-dimensional subspaces of
# F_q^{2d+1} which are not subspaces of a hyperplane H, and (d-1)-dimensional
# subspaces of H.
InstallMethod(TwistedGrassmannGraphCons, "as a spaces graph", true,
    [IsSpacesGraph, IsInt, IsInt], 0, function(filter, d, q)
        local F, H, S, V, n;
        n := 2*d+1;
        F := GF(q)^n;
        H := Subspace(F, CanonicalBasis(F){[1..n-1]}, "basis");
        V := Union(Filtered(Subspaces(F, d+1), x -> not IsSubset(H, x)),
                Subspaces(H, d-1));
        S := OnSubspaces(F);
        return Graph(Stabilizer(GL(n, q), H, S), V, S, function(x, y)
            local dx;
            dx := Dimension(x);
            if dx = Dimension(y) then
                return Dimension(Intersection(x, y)) = dx-1;
            else
                return IsSubset(x, y) or IsSubset(y, x);
            fi;
        end, true);
    end);

InstallMethod(TwistedGrassmannGraphCons, "default", true,
    [IsObject, IsInt, IsInt], 0, function(filter, d, q)
        return TwistedGrassmannGraphCons(IsSpacesGraph, d, q);
    end);

BindGlobal("TwistedGrassmannGraph", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j+1 then
        return TwistedGrassmannGraphCons(filt, arg[j], arg[j+1]);
    else
        Error("usage: TwistedGrassmannGraph( [<filter>, ]<int>, <int> )");
    fi;
end);

# The polar graph O^{(+/-)}(d, q) of isotropic lines of F_q^d
# with respect to a nondegenerate quadratic form.
InstallMethod(PolarGraphOCons, "as a spaces graph", true,
    [IsSpacesGraph, IsInt, IsInt, IsInt], 0, function(filter, e, d, q)
        local B, G, H, Q, V;
        G := GO(e, d, q);
        Q := InvariantQuadraticForm(G).matrix;
        B := Q + TransposedMat(Q);
        V := GF(q)^d;
        H := Graph(G, IsotropicSpacesQuadraticForm(Q)(Subspaces(V, 1)),
                OnSubspaces(V), function(x, y)
                    return x <> y and IsZero(Elements(x)[2]*B*Elements(y)[2]);
                end, true);
        H.duality := Sum;
        H.primality := Intersection;
        return H;
    end);

InstallMethod(PolarGraphOCons, "default", true,
    [IsObject, IsInt, IsInt, IsInt], 0, function(filter, e, d, q)
        return PolarGraphOCons(IsSpacesGraph, e, d, q);
    end);

BindGlobal("PolarGraphO", function(arg)
    local d, e, j, q, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j+1 then
        e := 0;
        d := arg[j];
        q := arg[j+1];
    elif Length(arg) = j+2 then
        e := arg[j];
        d := arg[j+1];
        q := arg[j+2];
    else
        Error("usage: PolarGraphO( [<filter>, ][<int>, ]<int>, <int> )");
    fi;
    return PolarGraphOCons(filt, e, d, q);
end);

# The polar graph NO^{+/-}orth(d, q) of nonisotropic lines of F_q^d
# with respect to a nondegenerate quadratic form.
InstallMethod(PolarGraphNOorthCons, "as a spaces graph", true,
    [IsSpacesGraph, IsInt, IsInt, IsInt], 0, function(filter, e, d, q)
        local f, z, B, G, H, Q, V;
        if d mod 2 = 1 then
            f := 0;
            if q mod 4 = 3 and d mod 4 = 1 then
                e := -e;
            fi;
        else
            f := e;
        fi;
        G := GO(f, d, q);
        Q := InvariantQuadraticForm(G).matrix;
        B := Q + TransposedMat(Q);
        V := GF(q)^d;
        z := Z(q)^((1-e)/2);
        H := Graph(G, NonisotropicSpacesQuadraticForm(Q, z)(Subspaces(V, 1)),
                OnSubspaces(V), function(x, y)
                    return x <> y and IsZero(Elements(x)[2]*B*Elements(y)[2]);
                end, true);
        return H;
    end);

InstallMethod(PolarGraphNOorthCons, "default", true,
    [IsObject, IsInt, IsInt, IsInt], 0, function(filter, e, d, q)
        return PolarGraphNOorthCons(IsSpacesGraph, e, d, q);
    end);

BindGlobal("PolarGraphNOorth", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j+2 then
        return PolarGraphNOorthCons(filt, arg[j], arg[j+1], arg[j+2]);
    else
        Error("usage: PolarGraphNOorth( [<filter>, ]<int>, <int>, <int> )");
    fi;
end);

# The polar graph Sp(d, q) of isotropic lines of F_q^d
# with respect to a nondegenerate symplectic form.
InstallMethod(PolarGraphSpCons, "as a spaces graph", true,
    [IsSpacesGraph, IsInt, IsInt], 0, function(filter, d, q)
        local G, H, Q, V;
        G := Sp(d, q);
        Q := InvariantBilinearForm(G).matrix;
        V := GF(q)^d;
        H := Graph(G, IsotropicSpacesQuadraticForm(Q)(Subspaces(V, 1)),
                OnSubspaces(V), function(x, y)
                    return x <> y and IsZero(Elements(x)[2]*Q*Elements(y)[2]);
                end, true);
        H.duality := Sum;
        H.primality := Intersection;
        return H;
    end);

InstallMethod(PolarGraphSpCons, "default", true,
    [IsObject, IsInt, IsInt], 0, function(filter, d, q)
        return PolarGraphSpCons(IsSpacesGraph, d, q);
    end);

BindGlobal("PolarGraphSp", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j+1 then
        return PolarGraphSpCons(filt, arg[j], arg[j+1]);
    else
        Error("usage: PolarGraphSp( [<filter>, ]<int>, <int> )");
    fi;
end);

# The polar graph U(d, r) of isotropic lines of F_{r^2}^d
# with respect to a nondegenerate Hermitean form.
InstallMethod(PolarGraphUCons, "as a spaces graph", true,
    [IsSpacesGraph, IsInt, IsInt], 0, function(filter, d, r)
        local c, B, F, G, H, Q, V;
        G := GU(d, r);
        Q := InvariantSesquilinearForm(G).matrix;
        V := GF(r^2)^d;
        B := Elements(CanonicalBasis(V));
        c := Conjugates(GF(r^2), GF(r), Z(r^2));
        F := x -> List(x, y -> y^r);
        H := Graph(G, IsotropicSpacesSesquilinearForm(Q, r)(Subspaces(V, 1)),
            OnSubspaces(V), function(x, y)
                return x <> y and IsZero(Elements(x)[2]*Q*F(Elements(y)[2]));
            end, true);
        H.duality := Sum;
        H.primality := Intersection;
        return H;
    end);

InstallMethod(PolarGraphUCons, "default", true,
    [IsObject, IsInt, IsInt], 0, function(filter, d, r)
        return PolarGraphUCons(IsSpacesGraph, d, r);
    end);

BindGlobal("PolarGraphU", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j+1 then
        return PolarGraphUCons(filt, arg[j], arg[j+1]);
    else
        Error("usage: PolarGraphU( [<filter>, ]<int>, <int> )");
    fi;
end);

# The dual polar graph B_d(q) of isotropic d-dimensional subspaces of
# F_q^{2d+1} with respect to a nondegenerate quadratic form.
InstallMethod(DualPolarGraphBCons,
    "as a spaces graph with full automorphism group", true,
    [IsSpacesGraph and FullAutomorphismGroup, IsInt, IsInt], 0,
    function(filter, d, q)
        local G, H, S, V, e, z, dp, invt;
        e := 2*d+1;
        V := GF(q)^e;
        G := GO(e, q);
        dp := DirectProduct(G, FieldExponentiationPermutationGroup(q));
        if d = 1 or q mod 2 = 1 then
            S := IsotropicSpacesQuadraticForm(InvariantQuadraticForm(G).matrix);
            if d = 1 then
                H := CompleteGraph(SymmetricGroup(q+1), q+1);
                AssignVertexNames(H, S(Subspaces(V, 1)));
                H.duality := Intersection;
                H.primality := Sum;
                return H;
            fi;
            invt := true;
        else
            S := [Subspace(V, Elements(CanonicalBasis(V)){[2..d+1]}, "basis")];
            invt := false;
        fi;
        return SubspaceGraph(dp, S, V, d, OnDualPolar(q, V, dp), invt);
    end);

InstallMethod(DualPolarGraphBCons, "as a spaces graph", true,
    [IsSpacesGraph, IsInt, IsInt], 0, function(filter, d, q)
        local G, S, V, e, z, invt;
        e := 2*d+1;
        V := GF(q)^e;
        G := GO(e, q);
        if q mod 2 = 1 then
            S := IsotropicSpacesQuadraticForm(InvariantQuadraticForm(G).matrix);
            invt := true;
        else
            S := [Subspace(V, Elements(CanonicalBasis(V)){[2..d+1]}, "basis")];
            invt := false;
        fi;
        return SubspaceGraph(G, S, V, d, invt);
    end);

InstallMethod(DualPolarGraphBCons, "with full automorphism group", true,
    [FullAutomorphismGroup, IsInt, IsInt], 0, function(filter, d, q)
        return DualPolarGraphBCons(IsSpacesGraph and FullAutomorphismGroup,
                                    d, q);
    end);

InstallMethod(DualPolarGraphBCons, "default", true,
    [IsObject, IsInt, IsInt], 0, function(filter, d, q)
        return DualPolarGraphBCons(IsSpacesGraph, d, q);
    end);

BindGlobal("DualPolarGraphB", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j+1 then
        return DualPolarGraphBCons(filt, arg[j], arg[j+1]);
    else
        Error("usage: DualPolarGraphB( [<filter>, ]<int>, <int> )");
    fi;
end);

# The dual polar graph C_d(q) of isotropic d-dimensional subspaces of
# F_q^{2d} with respect to a nondegenerate symplectic form.
InstallMethod(DualPolarGraphCCons,
    "as a spaces graph with full automorphism group", true,
    [IsSpacesGraph and FullAutomorphismGroup, IsInt, IsInt], 0,
    function(filter, d, q)
        local G, H, V, e, dp;
        e := 2*d;
        V := GF(q)^e;
        G := Sp(e, q);
        if d = 1 then
            H := CompleteGraph(SymmetricGroup(q+1), q+1);
            AssignVertexNames(H,
                IsotropicSpacesBilinearForm(InvariantBilinearForm(G).matrix)
                                            (Subspaces(V, 1)));
            H.duality := Intersection;
            H.primality := Sum;
            return H;
        fi;
        dp := DirectProduct(G, FieldExponentiationPermutationGroup(q));
        return SubspaceGraph(dp,
            [Subspace(V, Elements(CanonicalBasis(V)){[1..d]}, "basis")],
            V, d, OnDualPolar(q, V, dp), false);
    end);

InstallMethod(DualPolarGraphCCons, "as a spaces graph", true,
    [IsSpacesGraph, IsInt, IsInt], 0, function(filter, d, q)
        local V, e;
        e := 2*d;
        V := GF(q)^e;
        return SubspaceGraph(Sp(e, q),
            [Subspace(V, Elements(CanonicalBasis(V)){[1..d]}, "basis")],
            V, d, false);
    end);

InstallMethod(DualPolarGraphCCons, "with full automorphism group", true,
    [FullAutomorphismGroup, IsInt, IsInt], 0, function(filter, d, q)
        return DualPolarGraphCCons(IsSpacesGraph and FullAutomorphismGroup,
                                    d, q);
    end);

InstallMethod(DualPolarGraphCCons, "default", true,
    [IsObject, IsInt, IsInt], 0, function(filter, d, q)
        return DualPolarGraphCCons(IsSpacesGraph, d, q);
    end);

BindGlobal("DualPolarGraphC", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j+1 then
        return DualPolarGraphCCons(filt, arg[j], arg[j+1]);
    else
        Error("usage: DualPolarGraphC( [<filter>, ]<int>, <int> )");
    fi;
end);

# The dual polar graph C_d(q) of isotropic d-dimensional subspaces of
# F_q^{2d} with respect to a nondegenerate symplectic form.
InstallMethod(DualPolarGraphDCons,
    "as a spaces graph with full automorphism group", true,
    [IsSpacesGraph and FullAutomorphismGroup, IsInt, IsInt], 0,
    function(filter, d, q)
        local G, H, P, V, e, dp, invt;
        e := 2*d;
        V := GF(q)^e;
        G := GO(1, e, q);
        dp := DirectProduct(G, FieldExponentiationPermutationGroup(q));
        if d = 1 then
            return DualPolarGraphDCons(IsSpacesGraph, 1, q);
        elif d = 2 then
            H := CompleteMultipartiteGraph(2, q+1);
            P := IsotropicSpacesQuadraticForm(InvariantQuadraticForm(G).matrix)
                                                (Subspaces(V, 2));
            SortBy(P, v -> Dimension(Intersection(P[1], v)) mod 2);
            AssignVertexNames(H, P);
            H.duality := Intersection;
            H.primality := Sum;
            return H;
        elif q mod 2 = 0 then
            return SubspaceGraph(dp, [Subspace(V,
                            Elements(CanonicalBasis(V)){[1,3..e-1]}, "basis")],
                                    V, d, OnDualPolar(q, V, dp), false);
        else
            TryNextMethod();
        fi;
    end);

InstallMethod(DualPolarGraphDCons, "as a spaces graph", true,
    [IsSpacesGraph, IsInt, IsInt], 0, function(filter, d, q)
        local F, G, S, e, invt;
        e := 2*d;
        F := GF(q)^e;
        G := GO(1, e, q);
        if q mod 2 = 1 then
            S := IsotropicSpacesQuadraticForm(InvariantQuadraticForm(G).matrix);
            invt := true;
        else
            S := [Subspace(F, Elements(CanonicalBasis(F)){[1,3..e-1]},
                            "basis")];
            invt := false;
        fi;
        return SubspaceGraph(G, S, F, d, invt);
    end);

InstallMethod(DualPolarGraphDCons, "with full automorphism group", true,
    [FullAutomorphismGroup, IsInt, IsInt], 0, function(filter, d, q)
        return DualPolarGraphDCons(IsSpacesGraph and FullAutomorphismGroup,
                                    d, q);
    end);

InstallMethod(DualPolarGraphDCons, "default", true,
    [IsObject, IsInt, IsInt], 0, function(filter, d, q)
        return DualPolarGraphDCons(IsSpacesGraph, d, q);
    end);

# The dual polar graph D_d(q) of isotropic d-dimensional subspaces of
# F_q^{2d} with respect to a nondegenerate quadratic form of Witt index d.
BindGlobal("DualPolarGraphD", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j+1 then
        return DualPolarGraphDCons(filt, arg[j], arg[j+1]);
    else
        Error("usage: DualPolarGraphD( [<filter>, ]<int>, <int> )");
    fi;
end);

# The dual polar graph ^2D_{d+1}(q) of isotropic d-dimensional subspaces of
# F_q^{2d+2} with respect to a nondegenerate quadratic form of Witt index d.
InstallMethod(DualPolarGraph2DCons,
    "as a spaces graph with full automorphism group", true,
    [IsSpacesGraph and FullAutomorphismGroup, IsInt, IsInt], 0,
    function(filter, d, q)
        local G, H, S, V, e, dp;
        e := 2*d+2;
        V := GF(q)^e;
        G := GO(-1, e, q);
        dp := DirectProduct(G, FieldExponentiationPermutationGroup(q));
        S := IsotropicSpacesQuadraticForm(InvariantQuadraticForm(G).matrix);
        if d = 1 then
            H := CompleteGraph(SymmetricGroup(q^2+1), q^2+1);
            AssignVertexNames(H, S(Subspaces(V, 1)));
            H.duality := Intersection;
            H.primality := Sum;
            return H;
        elif q mod 2 = 0 and Log2Int(q) mod 2 = 1 then
            return SubspaceGraph(dp, S, V, d, OnDualPolar(q, V, dp), true);
        else
            TryNextMethod();
        fi;
    end);

InstallMethod(DualPolarGraph2DCons, "as a spaces graph", true,
    [IsSpacesGraph, IsInt, IsInt], 0, function(filter, d, q)
        local G, e;
        e := 2*d+2;
        G := GO(-1, e, q);
        return SubspaceGraph(G,
            IsotropicSpacesQuadraticForm(InvariantQuadraticForm(G).matrix),
                                GF(q)^e, d, true);
    end);

InstallMethod(DualPolarGraph2DCons, "with full automorphism group", true,
    [FullAutomorphismGroup, IsInt, IsInt], 0, function(filter, d, q)
        return DualPolarGraph2DCons(IsSpacesGraph and FullAutomorphismGroup,
                                    d, q);
    end);

InstallMethod(DualPolarGraph2DCons, "default", true,
    [IsObject, IsInt, IsInt], 0, function(filter, d, q)
        return DualPolarGraph2DCons(IsSpacesGraph, d, q);
    end);

BindGlobal("DualPolarGraph2D", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j+1 then
        return DualPolarGraph2DCons(filt, arg[j], arg[j+1]);
    else
        Error("usage: DualPolarGraph2D( [<filter>, ]<int>, <int> )");
    fi;
end);

# The dual polar graph ^2A_{e-1}(r) of isotropic [e/2]-dimensional
# subspaces of F_{r^2}^e with respect to a nondegenerate Hermitean form.
InstallMethod(DualPolarGraph2ACons,
    "as a spaces graph with full automorphism group", true,
    [IsSpacesGraph and FullAutomorphismGroup, IsInt, IsInt], 0,
    function(filter, e, r)
        local B, G, H, V, c, d, n, dp;
        V := GF(r^2)^e;
        B := Elements(CanonicalBasis(V));
        c := Conjugates(GF(r^2), GF(r), Z(r^2));
        d := Int(e/2);
        G := GU(e, r);
        dp := DirectProduct(G, FieldExponentiationPermutationGroup(r^2));
        if e in [2, 3] then
            n := r^(2*e-3) + 1;
            H := CompleteGraph(SymmetricGroup(n), n);
            AssignVertexNames(H,
        IsotropicSpacesSesquilinearForm(InvariantSesquilinearForm(G).matrix, r)
                                        (Subspaces(V, 1)));
            H.duality := Intersection;
            H.primality := Sum;
            return H;
        fi;
        return SubspaceGraph(dp,
            [Subspace(V, B{[1..d]} + (c[1]-c[2])*B{e+1-[1..d]}, "basis")],
                                V, d, OnDualPolar(r^2, V, dp), false);
    end);

InstallMethod(DualPolarGraph2ACons, "as a spaces graph", true,
    [IsSpacesGraph, IsInt, IsInt], 0, function(filter, e, r)
        local B, V, c, d;
        V := GF(r^2)^e;
        B := Elements(CanonicalBasis(V));
        c := Conjugates(GF(r^2), GF(r), Z(r^2));
        d := Int(e/2);
        return SubspaceGraph(GU(e, r),
            [Subspace(V, B{[1..d]} + (c[1]-c[2])*B{e+1-[1..d]}, "basis")],
                                V, d, false);
    end);

InstallMethod(DualPolarGraph2ACons, "with full automorphism group", true,
    [FullAutomorphismGroup, IsInt, IsInt], 0, function(filter, e, r)
        return DualPolarGraph2ACons(IsSpacesGraph and FullAutomorphismGroup,
                                    e, r);
    end);

InstallMethod(DualPolarGraph2ACons, "default", true,
    [IsObject, IsInt, IsInt], 0, function(filter, e, r)
        return DualPolarGraph2ACons(IsSpacesGraph, e, r);
    end);

BindGlobal("DualPolarGraph2A", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j+1 then
        return DualPolarGraph2ACons(filt, arg[j], arg[j+1]);
    else
        Error("usage: DualPolarGraph2A( [<filter>, ]<int>, <int> )");
    fi;
end);

# The Doro graph of nonisotropic 1-dimensional subspaces of F_q^4 with respect
# to a nondegenerate quadratic form. It is distance-regular for q = 4, 5.
InstallMethod(DoroGraphCons,
    "as a spaces graph with full automorphism group", true,
    [IsSpacesGraph and FullAutomorphismGroup, IsInt], 0, function(filter, q)
        local F, G, H, P, Q, V, dp;
        F := GF(q);
        V := F^4;
        G := GO(-1, 4, q);
        Q := InvariantQuadraticForm(G).matrix;
        dp := DirectProduct(G, FieldExponentiationPermutationGroup(q));
        if q <= 3 then
            H := NullGraph(SymmetricGroup(5*q));
            if q = 2 then
                AddEdgeOrbit(H, [1, 1]);
            fi;
            AssignVertexNames(H, NonisotropicSpacesQuadraticForm(Q,
                                                    Z(q)^0)(Subspaces(V, 1)));
            return H;
        elif LogInt(q, Characteristic(F)) mod 2 = 1 then
            return Graph(dp, [Subspace(V, BasisVectors(Basis(V)){[4]},
                            "basis")], OnDualPolar(q, V, dp), IsHyperbolic(Q));
        else
            TryNextMethod();
        fi;
    end);

InstallMethod(DoroGraphCons, "as a spaces graph", true,
    [IsSpacesGraph, IsInt], 0, function(filter, q)
        local G, V;
        V := GF(q)^4;
        G := GO(-1, 4, q);
        return Graph(G, [Subspace(V, BasisVectors(Basis(V)){[4]}, "basis")],
            OnSubspaces(V), IsHyperbolic(InvariantQuadraticForm(G).matrix));
    end);

InstallMethod(DoroGraphCons, "with full automorphism group", true,
    [FullAutomorphismGroup, IsInt], 0, function(filter, q)
        return DoroGraphCons(IsSpacesGraph and FullAutomorphismGroup, q);
    end);

InstallMethod(DoroGraphCons, "default", true, [IsObject, IsInt], 0,
    function(filter, q)
        return DoroGraphCons(IsSpacesGraph, q);
    end);

BindGlobal("DoroGraph", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j then
        return DoroGraphCons(filt, arg[j]);
    else
        Error("usage: DoroGraph( [<filter>, ]<int> )");
    fi;
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
