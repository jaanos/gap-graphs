# Complete multipartite graphs.
BindGlobal("CompleteMultipartiteGraph", function(arg)
    local sizes, dp, F, G;
    F := function(x, y) return x[1] <> y[1]; end;
    if Length(arg) = 1 then
        sizes := arg[1];
        dp := DirectProduct(List(sizes, i -> SymmetricGroup(i)));
        return Graph(dp, Union(List([1..Length(sizes)],
            i -> List([1..sizes[i]], j -> [i, j]))),
            OnSum(dp), F, true);
    else
        return ProductGraph([CompleteGraph(SymmetricGroup(arg[1])),
                            CompleteGraph(SymmetricGroup(arg[2]))], F);
    fi;
end);

# Cocktail party graphs.
BindGlobal("CocktailPartyGraph",
    n -> CompleteMultipartiteGraph(n, 2)
);

# Latin square graphs.
BindGlobal("LatinSquareGraph", function(G)
    local dim;
    if IsGroup(G) then
        return Graph(G, Cartesian(G, G), OnLatinSquare,
            function(x, y)
                return x <> y and (x[1] = y[1] or x[2] = y[2]
                                or x[1]*x[2] = y[1]*y[2]);
            end, true);
    else
        dim := DimensionsMat(G);
        return AdjFunGraph(Cartesian([1..dim[1]], [1..dim[2]]),
            function(x, y)
                return x <> y and (x[1] = y[1] or x[2] = y[2]
                                or G[x[1]][x[2]] = G[y[1]][y[2]]);
            end);
    fi;
end);
