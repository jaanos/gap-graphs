# The Hamming graph H(d, e) of vectors with d elements
# over an alphabet with e elements.
BindGlobal("HammingGraph", function(d, e)
    return Graph(WreathProduct(SymmetricGroup(e), SymmetricGroup(d)),
        Elements(ZmodnZ(e)^d), OnZmodnZVectors(d, e),
        function(x, y) return WeightVecFFE(x-y) = 1; end, true);
end);

# The Shrikhande graph with parameters v = 16, k = 6, lm = 2, mu = 2,
# i.e., the same as H(2, 4), but not isomorphic to it.
BindGlobal("ShrikhandeGraph",
    ComplementGraph(LatinSquareGraph(CyclicGroup(4))));
    
# The Doob graph Doob(n, d) of diameter 2*n+d
# as a box product of n copies of the Shrikhande graph and H(d, 4).
BindGlobal("DoobGraph", function(n, d)
    local l;
    l := ListWithIdenticalEntries(n, ShrikhandeGraph);
    if d > 0 then
        Add(l, HammingGraph(d, 4));
    fi;
    return BoxProductGraph(l);
end);
