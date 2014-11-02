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
    AssignVertexNames(G, Cartesian(List(Gs, H -> H.names)));
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
        local l;
        l := List([1..Length(Gs)], i -> Distance(Gs[i], x[i], y[i]));
        return WeightVecFFE(l) > 0 and Maximum(l) = 1;
    end);
end);

# The line graph of a graph.
BindGlobal("LineGraph",
    G -> Graph(G.group, Union(List([1..Length(G.representatives)],
        i -> List(G.adjacencies[i], y -> Set([G.representatives[i], y])))),
        OnSets, function(x, y) return Size(Intersection(x, y)) = 1; end)
);

# A graph with the set of d-dimensional subspaces P as the vertex set,
# with two subspaces being adjacent iff their intersection has dimension d-1.
BindGlobal("SubspaceGraph", function(P, d)
    return AdjFunGraph(P, function(x,y)
            return Dimension(Intersection(x,y)) = d-1;
        end);
end);

# The dual geometry graph of a geometric graph.
BindGlobal("DualGeometryGraph",
    G -> AdjFunGraph(Cliques(G), function(x,y)
            return Size(Intersection(x,y)) = 1;
        end)
);

# The flag graph of a point-line geometry
BindGlobal("FlagGraph",
    pg -> AdjFunGraph(Union(
            List([1..Length(pg)], x -> List(pg[x], y -> [x, y]))),
        function(x,y) return x <> y and (x[1] = y[1] or x[2] = y[2]); end)
);

# The incidence graph of a point-line geometry
BindGlobal("IncGraph",
    pg -> AdjFunGraph(Union(List([1..Length(pg)], x -> [0, x]),
            List(Union(pg), x -> [1, x])),
        function(x,y)
            local p, l;
            if x[1] = y[1] then
                return false;
            fi;
            if x[1] = 0 then
                l := x[2];
                p := y[2];
            else
                p := x[2];
                l := y[2];
            fi;
            return p in pg[l];
        end)
);

# The extended bipartite double of a graph.
BindGlobal("ExtendedBipartiteDoubleGraph",
    G -> Graph(G.group, Cartesian(G.representatives, GF(2)),
        function(x,y) return x; end,
        function(x,y)
            return x[2] <> y[2] and Distance(G, x[1], y[1]) <= 1;
        end, true)
);
