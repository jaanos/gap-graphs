LineGraph := function(G)
    return Graph(Group(()),
        Union(List([1..G.order],
            x -> List(Filtered(G.adjacencies[x], y -> y > x),
                z -> [x, z]))), function(x, y) return x; end,
        function(x, y) return Size(Intersection(x, y)) = 1; end, true);;
end;
