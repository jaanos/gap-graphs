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

# Cocktail party graphs
BindGlobal("CocktailPartyGraph",
    n -> CompleteMultipartiteGraph(n, 2)
);
