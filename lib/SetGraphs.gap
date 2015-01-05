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

# The folded Johnson graph.
BindGlobal("FoldedJohnsonGraph",
    d -> AntipodalQuotientGraph(JohnsonGraph(2*d, d))
);
