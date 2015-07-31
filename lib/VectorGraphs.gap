# The Hamming graph H(d, e) of vectors with d elements
# over an alphabet with e elements.
BindGlobal("HammingGraph", function(d, e)
    return Graph(WreathProductSymmetricGroups(e, d),
        Elements(ZmodnZ(e)^d), OnZmodnZVectors(d, e),
        function(x, y) return WeightVecFFE(x-y) = 1; end, true);
end);

# The d-dimensional hypercube
BindGlobal("HypercubeGraph", d -> HammingGraph(d, 2));

# The Shrikhande graph with parameters v = 16, k = 6, lm = 2, mu = 2,
# i.e., the same as H(2, 4), but not isomorphic to it.
BindGlobal("ShrikhandeGraph", ComplementGraph(LatinSquareGraph(Group(Z(5)))));
    
# The Doob graph Doob(n, d) of diameter 2*n+d
# as a box product of n copies of the Shrikhande graph and H(d, 4).
BindGlobal("DoobGraph", function(n, d)
    local l;
    l := ListWithIdenticalEntries(n, ShrikhandeGraph);
    if d > 0 then
        Add(l, HammingGraph(d, 4));
    fi;
    return BoxProductGraph(l);
end);

# The halved d-cube.
BindGlobal("HalvedCubeGraph",
    d -> HalvedGraph(HammingGraph(d, 2))
);

# The folded d-cube.
BindGlobal("FoldedCubeGraph",
    d -> AntipodalQuotientGraph(HammingGraph(d, 2))
);

# The folded halved 2d-cube.
BindGlobal("FoldedHalvedCubeGraph",
    d -> AntipodalQuotientGraph(HalvedGraph(HammingGraph(2*d, 2)))
);

# The Brouwer graph Br(q) of pairs of 3-dimensional vectors over F_q,
# with two pairs being adjacent whenever the difference of the first vectors
# equals the cross product of the second vectors.
BindGlobal("BrouwerGraph", function(q)
    return Graph(WreathProduct(FieldAdditionPermutationGroup(q),
            MatrixColumnEvenPermutationGroup(2, 3)), Elements(GF(q)^[2,3]),
        OnVectorPairs(q), function(x, y)
            return x <> y and x[1] - y[1] = VectorProduct(x[2], y[2]);
        end, true);
end);

# The Pasechnik graph Pa(q) as the extended bipartite double
# of the Brouwer graph Br(q).
BindGlobal("PasechnikGraph",
    q -> ExtendedBipartiteDoubleGraph(BrouwerGraph(q))
);

# The additive symplectic cover of the complete graph on q^{2n} vertices.
BindGlobal("AdditiveSymplecticCoverGraph", function(arg)
    local B, F, G, K, V, h, m, q, dp;
    if Length(arg) < 2 then
        Error("at least two arguments expected");
        return fail;
    fi;
    if Length(arg) > 2 then
        h := arg[3];
    else
        h := 0;
    fi;
    q := arg[1];
    F := GF(q);
    m := 2*arg[2];
    G := Sp(m, q);
    V := F^m;
    B := InvariantBilinearForm(G).matrix;
    dp := DirectProduct(Concatenation([G],
        ListWithIdenticalEntries(m+1, FieldAdditionPermutationGroup(q))));
    if h = 0 then
        K := AdditiveGroup(0*Z(q));
    else
        K := AdditiveGroup(BasisVectors(Basis(F)){[1..h]});
    fi;
    return Graph(dp, Cartesian(Unique(List(F, x -> x+K)), V),
        OnAdditiveSymplecticCover(q, m, B, dp),
        function(x, y)
            return x <> y and x[2]*B*y[2]+y[1]-Elements(x[1])[1] = K;
        end, true);
end);

# The multiplicative symplectic cover of the complete graph on q+1 vertices.
# It is distance-regular when m divides q-1 and either q or m is even.
BindGlobal("MultiplicativeSymplecticCoverGraph", function(q, m)
    local B, F, G, K, dp;
    F := GF(q);
    K := Group(Z(q)^((q-1)/m));
    G := Sp(2, q);
    B := InvariantBilinearForm(G).matrix;
    dp := DirectProduct(G, Group((1,2)),
                        FieldExponentiationPermutationGroup(q));
    return Graph(dp, Unique(List(Filtered(F^2, x -> x <> Zero(F^2)),
                                 v -> Set(List(K, g -> List(v, x -> g*x))))),
                 OnMultiplicativeSymplecticCover(q, dp), function(x, y)
                    return x <> y and x[1]*B*y[1] in K;
                 end, true);
end);

# The affine polar graph VO^{(+/-)}(d, q)
# with respect to a nondegenerate quadratic form.
BindGlobal("AffinePolarGraphVO", function(arg)
    local d, e, q, B, G, Q, V, dp;
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
    dp := DirectProduct(Concatenation([G],
            ListWithIdenticalEntries(d, FieldAdditionPermutationGroup(q))));
    Q := InvariantQuadraticForm(G).matrix;
    return Graph(dp, Elements(GF(q)^d), OnAffine(q, d, dp), function(x, y)
                    return x <> y and IsZero((x-y)*Q*(x-y));
                end, true);
end);

# The affine polar graph VNO^{(+/-)}(d, q)
# with respect to a nondegenerate quadratic form.
BindGlobal("AffinePolarGraphVNO", function(arg)
    local d, e, q, B, G, H, Q, V, dp;
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
    dp := DirectProduct(Concatenation([G],
            ListWithIdenticalEntries(d, FieldAdditionPermutationGroup(q))));
    Q := InvariantQuadraticForm(G).matrix;
    return Graph(dp, Elements(GF(q)^d), OnAffine(q, d, dp),
                function(x, y)
                    return x <> y and IsOne(((x-y)*Q*(x-y))^((q-1)/2));
                end, true);
end);
