# Complete multipartite graphs.
BindGlobal("CompleteMultipartiteGraph", function(arg)
    local sizes, dp, F, G, m, n;
    F := function(x, y) return x[1] <> y[1]; end;
    if Length(arg) = 0 then
        Error("at least one argument expected");
        return fail;
    elif Length(arg) = 1 then
        sizes := arg[1];
        dp := DirectProduct(List(sizes, i -> SymmetricGroup(i)));
        return Graph(dp, Union(List([1..Length(sizes)],
            i -> List([1..sizes[i]], j -> [i, j]))),
            OnSum(dp), F, true);
    else
        m := arg[1];
        n := arg[2];
        return Graph(WreathProductSymmetricGroups(n, m),
                     Cartesian([1..m], [1..n]),
                     function(x, g)
                        local y;
                        y := (n*(x[1]-1)+x[2])^g - 1;
                        return [Int(y/n)+1, y mod n + 1];
                     end, F, true);
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

# Hadamard graphs from Hadamard matrices.
BindGlobal("HadamardGraph", function(H)
    local F;
    if IsMatrix(H) then
        F := function(x)
            if x = H[1][1] then
                return 1;
            else
                return -1;
            fi;
        end;
        H := rec(matrix := List(H, l -> List(l, F)),
                 group := Group(()), transpose := ());
    fi;
    return Graph(Group(Union(GeneratorsOfGroup(H.group), [H.transpose])),
                 Cartesian([1, 2], [1, -1], [1..Length(H.matrix)]),
                 OnHadamardIndices(Length(H.matrix)),
                 function(x, y)
                    local p;
                    if x[1] = y[1] then
                        return false;
                    fi;
                    p := [];
                    p{[x[1], y[1]]} := [x[3], y[3]];
                    return H.matrix[p[1]][p[2]] = x[2]*y[2];
                 end, true);
end);
