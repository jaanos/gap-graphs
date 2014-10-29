RequirePackage("grape");

ExtendedBipartiteDouble := function(G)
    return Graph(G.group, Cartesian(G.representatives, GF(2)),
        function(x,y) return x; end,
        function(x,y)
            return x[2] <> y[2] and Distance(G, x[1], y[1]) <= 1;
        end, true);
end;
