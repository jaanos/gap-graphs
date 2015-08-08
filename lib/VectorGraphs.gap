# The Hamming graph H(d, e) of vectors with d elements
# over an alphabet with e elements.
InstallMethod(HammingGraphCons,
    "as a vector graph with full automorphism group", true,
    [IsVectorGraph and FullAutomorphismGroup, IsInt, IsInt], 0,
    function(filter, d, e)
        return BoxPowerGraph(NoVertexNames,
                             CompleteGraph(SymmetricGroup(e), e), d);
    end);

InstallMethod(HammingGraphCons, "as a vector graph", true,
    [IsVectorGraph, IsInt, IsInt], 0, function(filter, d, e)
        return HammingGraphCons(IsVectorGraph and FullAutomorphismGroup, d, e);
    end);

InstallMethod(HammingGraphCons, "with full automorphism group", true,
    [FullAutomorphismGroup, IsInt, IsInt], 0, function(filter, d, e)
        return HammingGraphCons(IsVectorGraph and FullAutomorphismGroup, d, e);
    end);

InstallMethod(HammingGraphCons, "default", true,
    [IsObject, IsInt, IsInt], 0, function(filter, d, e)
        return HammingGraphCons(IsVectorGraph, d, e);
    end);

BindGlobal("HammingGraph", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j+1 then
        return HammingGraphCons(filt, arg[j], arg[j+1]);
    else
        Error("usage: HammingGraph( [<filter>, ]<int>, <int> )");
    fi;
end);

# The d-dimensional hypercube
BindGlobal("HypercubeGraph", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j then
        return HammingGraphCons(filt, arg[j], 2);
    else
        Error("usage: HypercubeGraph( [<filter>, ]<int> )");
    fi;
end);
    
# The Doob graph Doob(n, d) of diameter 2*n+d
# as a box product of n copies of the Shrikhande graph and H(d, 4).
InstallMethod(DoobGraphCons,
    "as a vector graph with full automorphism group", true,
    [IsVectorGraph and FullAutomorphismGroup, IsInt, IsInt], 0,
    function(filter, n, d)
        return BoxProductGraph(
            BoxPowerGraph(ShrikhandeGraphCons(IsVectorGraph), n),
            HammingGraphCons(IsVectorGraph and FullAutomorphismGroup, d, 4));
    end);

InstallMethod(DoobGraphCons, "as a vector graph", true,
    [IsVectorGraph, IsInt, IsInt], 0, function(filter, n, d)
        return DoobGraphCons(IsVectorGraph and FullAutomorphismGroup, n, d);
    end);

InstallMethod(DoobGraphCons, "with full automorphism group", true,
    [FullAutomorphismGroup, IsInt, IsInt], 0, function(filter, n, d)
        return DoobGraphCons(IsVectorGraph and FullAutomorphismGroup, n, d);
    end);

InstallMethod(DoobGraphCons, "default", true,
    [IsObject, IsInt, IsInt], 0, function(filter, n, d)
        return DoobGraphCons(IsVectorGraph, n, d);
    end);

BindGlobal("DoobGraph", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j+1 then
        return DoobGraphCons(filt, arg[j], arg[j+1]);
    else
        Error("usage: DoobGraph( [<filter>, ]<int>, <int> )");
    fi;
end);

# The halved d-cube.
BindGlobal("HalvedCubeGraph", function(arg)
    local G, j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j then
        if arg[j] = 4 then
            G := CocktailPartyGraph(4);
            AssignVertexNames(G, List(G.names,
                t -> List([1..4], function(i)
                                    if t[1] = 1 or i = 1 or i = t[1] then
                                        return t[2];
                                    else
                                        return 3 - t[2];
                                    fi;
                                  end)));
            return G;
        else
            return HalvedGraph(HammingGraphCons(filt, arg[j], 2));
        fi;
    else
        Error("usage: HalvedCubeGraph( [<filter>, ]<int> )");
    fi;
end);

# The folded d-cube.
BindGlobal("FoldedCubeGraph", function(arg)
    local G, j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j then
        if arg[j] = 4 then
            G := CompleteMultipartiteGraph(2, 4);
            AssignVertexNames(G, List(List(G.names, t -> List([1..4],
                        function(i)
                            if t[1] = 1 and (t[2] = 1 or i = 1 or i = t[2])
                                or t[1] = 2 and i <> t[2] then
                                    return 1;
                            else
                                    return 2;
                            fi;
                        end)), s -> Set([s, 3-s])));
            return G;
        else
            return AntipodalQuotientGraph(HammingGraphCons(filt, arg[j], 2));
        fi;
    else
        Error("usage: FoldedCubeGraph( [<filter>, ]<int> )");
    fi;
end);

# The folded halved 2d-cube.
BindGlobal("FoldedHalvedCubeGraph", function(arg)
    local G, j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j then
        if arg[j] = 3 then
            G := CompleteGraph(SymmetricGroup(16));
            AssignVertexNames(G, List(Tuples([1,2], 4), function(t)
                    if WeightVecFFE(t-1) mod 2 = 0 then
                        return [Concatenation([1,1], t),
                                Concatenation([2,2], 3-t)];
                    else
                        return [Concatenation([1,2], t),
                                Concatenation([2,1], 3-t)];
                    fi;
                end));
            return G;
        else
            return AntipodalQuotientGraph(HalvedGraph(HammingGraphCons(filt,
                                                                2*arg[j], 2)));
        fi;
    else
        Error("usage: FoldedHalvedCubeGraph( [<filter>, ]<int> )");
    fi;
end);

DeclareSynonym("HalvedFoldedCubeGraph", FoldedHalvedCubeGraph);

# The Brouwer graph Br(q) of pairs of 3-dimensional vectors over F_q,
# with two pairs being adjacent whenever the difference of the first vectors
# equals the cross product of the second vectors.
InstallMethod(BrouwerGraphCons,
    "as a vector graph with full automorphism group", true,
    [IsVectorGraph and FullAutomorphismGroup, IsInt], 0, function(filter, q)
        local G, dp, wp;
        if q = 2 then
            G := FoldedCubeGraph(IsVectorGraph and FullAutomorphismGroup, 7);
            AssignVertexNames(G, List(G.names, function(p)
                                    local t, u, v;
                                    v := First(p, s -> WeightVecFFE(s-1) <= 3);
                                    t := List(Filtered([1..7], i -> v[i] = 2),
                                                j -> IntToBits(j, 3)*Z(2)^0);
                                    if Length(t) <= 1 then
                                        u := [0*Z(2), 0*Z(2), 0*Z(2)];
                                    elif Length(t) = 2 then
                                        u := VectorProduct(t[1], t[2]);
                                    else
                                        u := VectorProduct(t[1], t[2])
                                            + VectorProduct(t[1], t[3])
                                            + VectorProduct(t[2], t[3]);
                                    fi;
                                    if Length(t) = 0 then
                                        t[1] := u;
                                    fi;
                                    return [u, Sum(t)];
                                end));
            return G;
        else
            wp := WreathProduct(FieldAdditionPermutationGroup(q),
                                MatrixColumnEvenPermutationGroup(2, 3));
            dp := DirectProduct(wp, GL(3, q),
                                FieldExponentiationPermutationGroup(q));
            return Graph(dp, Elements(GF(q)^[2,3]), OnVectorPairs(q, dp, wp),
                function(x, y)
                    return x <> y and x[1] - y[1] = VectorProduct(x[2], y[2]);
                end, true);
        fi;
    end);

InstallMethod(BrouwerGraphCons, "as a vector graph", true,
    [IsVectorGraph, IsInt], 0, function(filter, q)
        local dp, wp;
        wp := WreathProduct(FieldAdditionPermutationGroup(q),
                            MatrixColumnEvenPermutationGroup(2, 3));
        dp := DirectProduct(wp, Group(IdentityMat(3, GF(q))), Group(()));
        return Graph(dp, Elements(GF(q)^[2,3]), OnVectorPairs(q, dp, wp),
            function(x, y)
                return x <> y and x[1] - y[1] = VectorProduct(x[2], y[2]);
            end, true);
    end);

InstallMethod(BrouwerGraphCons, "with full automorphism group", true,
    [FullAutomorphismGroup, IsInt], 0, function(filter, q)
        return BrouwerGraphCons(IsVectorGraph and FullAutomorphismGroup, q);
    end);

InstallMethod(BrouwerGraphCons, "as a vector graph", true,
    [IsObject, IsInt], 0, function(filter, q)
        return BrouwerGraphCons(IsVectorGraph, q);
    end);

BindGlobal("BrouwerGraph", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j then
        return BrouwerGraphCons(filt, arg[j]);
    else
        Error("usage: BrouwerGraph( [<filter>, ]<int> )");
    fi;
end);

# The Pasechnik graph Pa(q) as the extended bipartite double
# of the Brouwer graph Br(q).
InstallMethod(PasechnikGraphCons, "as a vector graph", true,
    [IsVectorGraph, IsInt], 0, function(filter, q)
        return ExtendedBipartiteDoubleGraph(BrouwerGraphCons(IsVectorGraph,
                                                             q));
    end);

InstallMethod(PasechnikGraphCons, "as a vector graph", true,
    [IsObject, IsInt], 0, function(filter, q)
        return PasechnikGraphCons(IsVectorGraph, q);
    end);

BindGlobal("PasechnikGraph", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j then
        return PasechnikGraphCons(filt, arg[j]);
    else
        Error("usage: PasechnikGraph( [<filter>, ]<int> )");
    fi;
end);

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
