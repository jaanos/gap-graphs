# The graph obtained from an adjacency function on the vertex set.
BindGlobal("AdjFunGraph", function(E, F)
    return Graph(Group(()), E, function(x, y) return x; end, F, true);
end);

# A generic product graph.
BindGlobal("ProductGraph", function(Gs, F)
    local G, GG, dp;
    dp := DirectProduct(List(Gs, H -> H.group));
    G := Graph(dp, Cartesian(List(Gs, H -> [1..H.order])),
        OnProduct(Length(Gs), dp), F, true);
    for GG in Gs do
        if not "names" in RecNames(GG) then
            GG.names := [1..GG.order];
        fi;
    od;
    AssignVertexNames(G, List(G.names,
        f -> List([1..Length(f)], i -> Gs[i].names[f[i]])));
    return G;
end);

# The box product of two or more graphs.
BindGlobal("BoxProductGraph", function(arg)
    local Gs;
    if Length(arg) = 1 then
        Gs := arg[1];
    else
        Gs := arg;
    fi;
    return ProductGraph(Gs, function(x, y)
        local l;
        l := List([1..Length(Gs)], i -> Distance(Gs[i], x[i], y[i]));
        return WeightVecFFE(l) = 1 and Sum(l) = 1;
    end);
end);

# The cross product of two or more graphs.
BindGlobal("CrossProductGraph", function(arg)
    local Gs;
    if Length(arg) = 1 then
        Gs := arg[1];
    else
        Gs := arg;
    fi;
    return ProductGraph(Gs, function(x, y)
        local l;
        l := List([1..Length(Gs)], i -> Distance(Gs[i], x[i], y[i]));
        return Minimum(l) = 1 and Maximum(l) = 1;
    end);
end);

# The strong product of two or more graphs.
BindGlobal("StrongProductGraph", function(arg)
    local Gs;
    if Length(arg) = 1 then
        Gs := arg[1];
    else
        Gs := arg;
    fi;
    return ProductGraph(Gs, function(x, y)
        return Maximum(List([1..Length(Gs)],
            i -> Distance(Gs[i], x[i], y[i]))) = 1;
    end);
end);

# The bipartite double of a graph.
BindGlobal("BipartiteDoubleGraph", function(G)
    local H;
    H := BipartiteDouble(G);
    CheckDualityFunctions(G);
    H.halfDuality := BipartiteDoubleDualityFunction(G.duality);
    H.halfPrimality := BipartiteDoubleDualityFunction(G.primality);
    return H;
end);

# The extended bipartite double of a graph.
BindGlobal("ExtendedBipartiteDoubleGraph", function(G)
    local dp, signs, H;
    signs := ["+", "-"];
    dp := DirectProduct(G.group, SymmetricGroup(2));
    H := Graph(dp, Cartesian(G.representatives, signs),
        OnSignedPoints(dp, signs), function(x, y)
            return x[2] <> y[2] and Distance(G, x[1], y[1]) <= 1;
        end);
    if "names" in RecNames(G) then
        AssignVertexNames(H, List(H.names, x -> [G.names[x[1]], x[2]]));
    fi;
    CheckDualityFunctions(G);
    H.halfDuality := BipartiteDoubleDualityFunction(G.duality);
    H.halfPrimality := BipartiteDoubleDualityFunction(G.primality);
    return H;
end);

# The halved graph of a bipartite graph. The optional second argument
# allows choosing between the first and second halves.
BindGlobal("HalvedGraph", function(arg)
    local n, G, G2, H, vs;
    if Length(arg) < 1 then
        Error("at least one argument expected");
        return fail;
    fi;
    G := arg[1];
    if Length(arg) > 1 then
        n := arg[2];
    else
        n := 1;
    fi;
    if not IsConnectedGraph(G) then
        Error("not a connected graph");
        return fail;
    fi;
    if not IsBipartite(G) then
        Error("not a bipartite graph");
        return fail;
    fi;
    G2 := DistanceGraph(G, 2);
    vs := ConnectedComponent(G2, 1);
    if n = 2 then
        vs := Difference([1..G.order], vs);
    fi;
    H := Graph(Stabilizer(G.group, vs, OnSets), vs, OnPoints,
        function(x, y) return Distance(G2, x, y) = 1; end, true);
    if "names" in RecNames(G) then
        AssignVertexNames(H, G.names{vs});
    fi;
    if "halfDuality" in RecNames(G) then
        if n = 2 then
            H.primality := G.halfDuality;
        else
            H.duality := G.halfDuality;
        fi;
    fi;
    if "halfPrimality" in RecNames(G) then
        if n = 2 then
            H.duality := G.halfPrimality;
        else
            H.primality := G.halfPrimality;
        fi;
    fi;
    return H;
end);

# The antipodal quotient of an antipodal cover.
BindGlobal("AntipodalQuotientGraph", function(G)
    local d, H;
    if not IsAntipodalCover(G) then
        Error("not an antipodal cover");
        return fail;
    fi;
    d := Diameter(G);
    H := Graph(G.group,
        Set(List(G.representatives, x -> DistanceSet(G, [0, d], x))),
        OnSets, function(x, y)
            return IsSubset(DistanceSet(G, 1, x), y);
        end);
    if "names" in RecNames(G) then
        AssignVertexNames(H, List(H.names, f -> G.names{f}));
    fi;
    return H;
end);

# A graph with the set of d-dimensional subspaces of V filtered by S
# as the vertex set, acted upon by the matrix group G,
# with two subspaces being adjacent iff their intersection has dimension d-1.
BindGlobal("SubspaceGraph", function(arg)
    local G, H, S, V, d, invt, vcs;
    if Length(arg) < 4 then
        Error("at least four arguments expected");
        return fail;
    fi;
    G := arg[1];
    S := arg[2];
    V := arg[3];
    d := arg[4];
    if Length(arg) > 4 then
        invt := arg[5];
    else
        invt := true;
    fi;
    if IsList(S) then
        vcs := S;
    else
        vcs := S(Subspaces(V, d));
    fi;
    H := Graph(G, vcs, OnSubspaces(V), function(x,y)
                    return Dimension(Intersection(x,y)) = d-1;
                end, invt);
    H.duality := Intersection;
    H.primality := Sum;
    return H;
end);

# The clique (dual geometry) graph of a collinearity graph. The optional second
# argument allows choosing a connected component of the resulting graph.
BindGlobal("CliqueGraph", function(arg)
    local C, G, H, n;
    if Length(arg) < 1 then
        Error("at least one argument expected");
        return fail;
    fi;
    G := arg[1];
    C := Cliques(G);
    if Length(arg) > 1 then
        n := arg[2];
        if not IsList(n) then
            n := [n];
        fi;
    else
        n := [1..Length(C)];
    fi;
    H := Graph(G.group, C{n}, OnSets,
                function(x, y)
                    return Size(Intersection(x,y)) = 1;
                end);
    if "names" in RecNames(G) then
        CheckDualityFunctions(G);
        H.duality := G.primality;
        H.primality := G.duality;
        AssignVertexNames(H, List(H.names, f -> G.duality(G.names{f})));
        if "halfDuality" in RecNames(G) then
            H.halfDuality := G.halfDuality;
        fi;
        if "halfPrimality" in RecNames(G) then
            H.halfPrimality := G.halfPrimality;
        fi;
    fi;
    return H;
end);

# The incidence graph of a collinearity graph.
BindGlobal("IncidenceGraph", function(arg)
    local C, G, H, n;
    if Length(arg) < 1 then
        Error("at least one argument expected");
        return fail;
    fi;
    G := arg[1];
    C := Cliques(G);
    if Length(arg) > 1 then
        n := arg[2];
        if not IsList(n) then
            n := [n];
        fi;
    else
        n := [1..Length(C)];
    fi;
    H := Graph(G.group, Union(G.representatives, C{n}),
                OnPointsOrLines(OnPoints, IsList),
                function(x, y)
                    return (IsList(y) and x in y) or (IsList(x) and y in x);
                end);
    if "names" in RecNames(G) and G.order > 0 then
        CheckDualityFunctions(G);
        if IsList(H.names[1]) then
            H.halfDuality := G.primality;
            H.halfPrimality := G.duality;
        else
            H.halfDuality := G.duality;
            H.halfPrimality := G.primality;
        fi;
        AssignVertexNames(H, List(H.names, function(f)
                                            if IsList(f) then
                                                return G.duality(G.names{f});
                                            else
                                                return G.names[f];
                                            fi;
        end));
    fi;
    return H;
end);
