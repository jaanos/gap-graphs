# The Kneser graph on k-subsets of a set with n elements.
BindGlobal("KneserGraph", function(arg)
    local G, n, k, invt, vcs;
    if Length(arg) < 2 then
        Error("at least two arguments expected");
        return fail;
    else
        n := arg[1];
        k := arg[2];
        if Length(arg) > 2 and 2*k < n then
            invt := arg[3];
        else
            invt := true;
        fi;
    fi;
    if n = 2*k then
        vcs := List(Combinations([2..n], k),
                    s -> [Difference([1..n], s), s]);
        G := ComplementGraph(CocktailPartyGraph(Length(vcs)));
        AssignVertexNames(G, List(G.names, t -> vcs[t[1]][t[2]]));
    else
        if invt then
            vcs := Combinations([1..n], k);
        else
            vcs := [[1..k]];
        fi;
        if n < 2*k then
            G := NullGraph(SymmetricGroup(Length(vcs)), Length(vcs));
            AssignVertexNames(G, vcs);
        else
            G := Graph(SymmetricGroup(n), vcs, OnSets, DisjointSets, invt);
        fi;
    fi;
    return G;
end);

# The Odd graph of diameter d on 2*d+1 points.
BindGlobal("OddGraph", d -> KneserGraph(2*d+1, d, false));

# The doubled Odd graph on 2*d+1 points.
BindGlobal("DoubledOddGraph", function(d)
    local n, dp;
    n := 2*d+1;
    dp := DirectProduct(SymmetricGroup(n), SymmetricGroup(2));
    return Graph(dp, Union(Combinations([1..n], d), Combinations([1..n], d+1)),
                    OnDoubledOdd(n, dp), SymmetrizedInclusion, true);
end);

# The Johnson graph on d-subsets of a set with n elements.
InstallMethod(JohnsonGraphCons, "as a set graph with full automorphism group",
    true, [IsSetGraph and FullAutomorphismGroup, IsInt, IsInt], 0,
    function(filter, n, d)
        local dp;
        if n = 2*d then
            dp := DirectProduct(SymmetricGroup(n), Group((1,2)));
            return Graph(dp, Combinations([1..n], d), OnJohnson(n, dp),
                         SetIntersection(d-1), true);
        else
            return JohnsonGraphCons(IsSetGraph, n, d);
        fi;
    end);

InstallMethod(JohnsonGraphCons, "as a set graph", true,
    [IsSetGraph, IsInt, IsInt], 0, function(filter, n, d)
        return JohnsonGraph(n, d);
    end);

InstallMethod(JohnsonGraphCons, "with full automorphism group", true,
    [FullAutomorphismGroup, IsInt, IsInt], 0, function(filter, n, d)
        return JohnsonGraphCons(IsSetGraph and FullAutomorphismGroup, n, d);
    end);

InstallMethod(JohnsonGraphCons, "default", true,
    [IsObject, IsInt, IsInt], 0, function(filter, n, d)
        return JohnsonGraphCons(IsSetGraph, n, d);
    end);

# The folded Johnson graph.
InstallMethod(FoldedJohnsonGraphCons,
    "as a set graph with full automorphism group", true,
    [IsSetGraph and FullAutomorphismGroup, IsInt], 0, function(filter, d)
        local G;
        if d = 3 then
            G := CompleteGraph(SymmetricGroup(10));
            AssignVertexNames(G, List(Combinations([2..6], 3),
                                        s -> [Difference([1..6], s), s]));
            return G;
        else
            return FoldedJohnsonGraphCons(IsSetGraph, d);
        fi;
    end);

InstallMethod(FoldedJohnsonGraphCons, "as a set graph", true,
    [IsSetGraph, IsInt], 0, function(filter, d)
        return AntipodalQuotientGraph(JohnsonGraphCons(IsSetGraph, 2*d, d));
    end);

InstallMethod(FoldedJohnsonGraphCons, "with full automorphism group", true,
    [FullAutomorphismGroup, IsInt], 0, function(filter, d)
        return FoldedJohnsonGraphCons(IsSetGraph and FullAutomorphismGroup, d);
    end);

InstallMethod(FoldedJohnsonGraphCons, "default", true,
    [IsObject, IsInt], 0, function(filter, d)
        return FoldedJohnsonGraphCons(IsSetGraph, d);
    end);

BindGlobal("FoldedJohnsonGraph", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j then
        return FoldedJohnsonGraphCons(filt, arg[j]);
    else
        Error("usage: FoldedJohnsonGraph( [<filter>, ]<int> )");
    fi;
end);

# The three Chang graphs with v=28, k=12, lm=6, mu=4
InstallMethod(ChangGraphCons,
    "as a set graph with full automorphism group", true,
    [IsSetGraph and FullAutomorphismGroup, IsInt], 0, function(filter, j)
        local J, dp;
        J := JohnsonGraphCons(IsSetGraph, 8, 2);
        dp := DirectProduct(Stabilizer(SymmetricGroup(8),
                                    ChangGraphSwitchingSet[j], OnSetsSets),
                            Group(ChangGraphInvolution[j]));
        return SwitchedGraph(J, List(ChangGraphSwitchingSet[j],
                                        x -> Position(J.names, x)),
                                Action(dp, J.names, OnChang(dp)));
    end);

InstallMethod(ChangGraphCons, "as a set graph", true,
    [IsSetGraph, IsInt], 0, function(filter, j)
        local J, S;
        J := JohnsonGraphCons(IsSetGraph, 8, 2);
        S := List(ChangGraphSwitchingSet[j], x -> Position(J.names, x));
        return SwitchedGraph(J, S, Stabilizer(J.group, S, OnSets));
    end);

InstallMethod(ChangGraphCons, "with full automorphism group", true,
    [FullAutomorphismGroup, IsInt], 0, function(filter, j)
        return ChangGraphCons(IsSetGraph and FullAutomorphismGroup, j);
    end);

InstallMethod(ChangGraphCons, "default", true,
    [IsObject, IsInt], 0, function(filter, j)
        return ChangGraphCons(IsSetGraph, j);
    end);

BindGlobal("ChangGraph", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j then
        return ChangGraphCons(filt, arg[j]);
    else
        Error("usage: ChangGraph( [<filter>, ]<int> )");
    fi;
end);
