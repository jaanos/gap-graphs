# The Kneser graph on k-subsets of a set with n elements.
BindGlobal("KneserGraph", function(arg)
    local n, k, invt, vcs;
    if Length(arg) < 2 then
        Error("at least two arguments expected");
        return fail;
    else
        n := arg[1];
        k := arg[2];
        if Length(arg) > 2 then
            invt := arg[3];
        else
            invt := true;
        fi;
    fi;
    if invt then
        vcs := Combinations([1..n], k);
    else
        vcs := [[1..k]];
    fi;
    return Graph(SymmetricGroup(n), vcs, OnSets, DisjointSets);
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

# The folded Johnson graph.
BindGlobal("FoldedJohnsonGraph",
    d -> AntipodalQuotientGraph(JohnsonGraph(2*d, d))
);

# The three Chang graphs with v=28, k=12, lm=6, mu=4
BindGlobal("ChangGraph", function(j)
    local J, S, Ss;
    Ss := [List([1..4], i -> [i, i+4]),
           Set(List([1..8], i -> Set([i, (i mod 8)+1]))),
           Union(List([1..3], i -> Set([i, (i mod 3)+1])),
                 List([1..5], i -> Set([i+3, (i mod 5)+4])))];
    J := JohnsonGraph(8, 2);
    S := List(Ss[j], x -> Position(J.names, x));
    return SwitchedGraph(J, S, Stabilizer(J.group, S, OnSets));
end);
