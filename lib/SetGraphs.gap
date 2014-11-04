# The Kneser graph on k-subsets of a set with n elements.
BindGlobal("KneserGraph", function(n, k)
    return Graph(SymmetricGroup(n), [[1..k]], OnSets,
        function(x,y) return Intersection(x,y)=[]; end);
end);

# The Odd graph of diameter d on 2*d+1 points.
BindGlobal("OddGraph", d -> KneserGraph(2*d+1, d));

# The folded Johnson graph.
BindGlobal("FoldedJohnsonGraph",
    d -> AntipodalQuotientGraph(JohnsonGraph(2*d, d))
);
