# The Hamming graph H(d, e) of vectors with d elements
# over an alphabet with e elements.
BindGlobal("HammingGraph", function(d, e)
    return Graph(WreathProduct(SymmetricGroup(e), SymmetricGroup(d)),
        Elements(ZmodnZ(e)^d), OnZmodnZVectors(d, e),
        function(x, y) return WeightVecFFE(x-y) = 1; end, true);
end);
