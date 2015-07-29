# Complete multipartite graphs given a list of part sizes
# with the full automorphism group.
InstallOtherMethod(CompleteMultipartiteGraphCons,
    "for a list of part sizes with full automorphism group", true,
    [FullAutomorphismGroup, IsList], 0, function(filter, sizes)
        local H, Gs, dp, l, p;
        l := Set(sizes);
        p := List(l, i -> Filtered([1..Length(sizes)], j -> sizes[j] = i));
        Gs := List([1..Length(l)],
                    i -> CompleteMultipartiteGraphCons(FullAutomorphismGroup,
                                                        Length(p[i]), l[i]));
        H := GraphJoin(NoVertexNames, Gs);
        AssignVertexNames(H, List(H.names, function(t)
                                                local tt;
                                                tt := Gs[t[1]].names[t[2]];
                                                return [p[t[1]][tt[1]], tt[2]];
                                            end));
        return H;
    end);

# Complete multipartite graphs given a list of part sizes,
InstallOtherMethod(CompleteMultipartiteGraphCons, "for a list of part sizes",
    true, [IsObject, IsList], 0, function(filter, sizes)
        local dp;
        dp := DirectProduct(List(sizes, i -> SymmetricGroup(i)));
        return Graph(dp, Union(List([1..Length(sizes)],
            i -> List([1..sizes[i]], j -> [i, j]))),
            OnSum(dp), DifferentParts, true);
    end);

InstallMethod(CompleteMultipartiteGraphCons,
    "for equal part sizes with full automorphism group", true,
    [FullAutomorphismGroup, IsInt, IsInt], 0, function(filter, m, n)
        return Graph(WreathProductSymmetricGroups(n, m),
                     Cartesian([1..m], [1..n]),
                     function(x, g)
                        local y;
                        y := (n*(x[1]-1)+x[2])^g - 1;
                        return [Int(y/n)+1, y mod n + 1];
                     end, DifferentParts, true);
    end);

InstallMethod(CompleteMultipartiteGraphCons, "for equal part sizes",
    true, [IsObject, IsInt, IsInt], 0, function(filter, m, n)
        return CompleteMultipartiteGraphCons(FullAutomorphismGroup, m, n);
    end);

BindGlobal("CompleteMultipartiteGraph", function(arg)
    local j, filt;
    if IsFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j then
        return CompleteMultipartiteGraphCons(filt, arg[j]);
    elif Length(arg) = j+1 then
        return CompleteMultipartiteGraphCons(filt, arg[j], arg[j+1]);
    else
        Error("usage: CompleteMultipartiteGraph( [<filter>, ]{<list> |<int>, <int> })");
    fi;
end);

# Cycle graphs.
BindGlobal("CycleGraph", n -> Graph(DihedralGroup(IsPermGroup, 2*n), [1..n],
    OnPoints, function(x, y)
        return (x-y) mod n in [1,n-1];
    end, true));

# Cocktail party graphs.
BindGlobal("CocktailPartyGraph",
    n -> CompleteMultipartiteGraph(n, 2)
);

# Paley graphs.
# For q = 1 (mod 4) a prime power, the graph is strongly regular.
# For q = 3 (mod 4) a prime power, the graph is directed.
BindGlobal("PaleyGraph", function(q)
    local dp;
    dp := DirectProduct(FieldAdditionPermutationGroup(q),
        Group(GeneratorsOfGroup(FieldMultiplicationPermutationGroup(q))[1]^2));
    return Graph(dp, Elements(GF(q)), OnPaley(q, dp),
        function(x, y)
            return IsOne((x-y)^((q-1)/2));
        end, true);
end);

# Latin square graphs.
BindGlobal("LatinSquareGraph", function(arg)
    local dim, dp, G, invt, vcs;
    if Length(arg) = 0 then
        Error("at least one argument expected");
        return fail;
    fi;
    G := arg[1];
    if Length(arg) > 1 then
        invt := arg[2];
    else
        invt := true;
    fi;
    if IsGroup(G) then
        dp := DirectProduct(G, G);
        if invt then
            vcs := Cartesian(G, G);
        else
            vcs := [[One(G), One(G)]];
        fi;
        return Graph(dp, vcs, OnLatinSquare(dp),
            function(x, y)
                return x <> y and (x[1] = y[1] or x[2] = y[2]
                                or x[1]*x[2] = y[1]*y[2]);
            end, invt);
    else
        dim := DimensionsMat(G);
        return AdjFunGraph(Cartesian([1..dim[1]], [1..dim[2]]),
            function(x, y)
                return x <> y and (x[1] = y[1] or x[2] = y[2]
                                or G[x[1]][x[2]] = G[y[1]][y[2]]);
            end);
    fi;
end);

# Complete Taylor graphs, i.e. complete bipartite graphs minus a matching.
BindGlobal("CompleteTaylorGraph", function(n)
    local G;
    G := EdgeOrbitsGraph(Group([(1,2)(n+1,n+2),
                    PermList(Concatenation([2..n], [1], [n+2..2*n], [n+1])),
                    PermList(Concatenation([n+1..2*n], [1..n]))]), [1, n+2]);
    AssignVertexNames(G, Cartesian([1, 2], [1..n]));
    return G;
end);
