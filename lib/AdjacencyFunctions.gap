# Adjacency function for complete multipartite graphs
BindGlobal("DifferentParts", function(x, y)
    return x[1] <> y[1];
end);

# Adjacency function graph joins
BindGlobal("GraphJoinAdjacency", Gs -> function(x, y)
    return x[1] <> y[1] or (x[1] = y[1] and
                            IsVertexPairEdge(Gs[x[1]], x[2], y[2]));
end);

# Adjacency function for Kneser-type graphs.
BindGlobal("DisjointSets", function(x, y)
    return Intersection(x, y) = [];
end);

# Adjacency function for conjugacy class graphs.
BindGlobal("GroupIntersection", n -> function(x, y)
    return Order(Intersection(x, y)) = n;
end);

# Adjacency function for set graphs.
BindGlobal("SetIntersection", n -> function(x, y)
    return Length(Intersection(x, y)) = n;
end);

# Adjacency function for doubled Odd and Grassmann graphs.
BindGlobal("SymmetrizedInclusion", function(x, y)
    return x <> y and (IsSubset(x, y) or IsSubset(y, x));
end);

# Adjacency function for roots of E_8.
BindGlobal("RootAdjacency", function(x, y)
    return x*y = 8;
end);

# Adjacency function for forms graphs.
BindGlobal("RankAdjacency",
    r -> function(x, y)
        return RankMat(x-y) in r;
    end);

# Adjacency function for adjacency matrices.
BindGlobal("MatrixAdjacency",
    M -> function(x, y)
        return M[x][y] = 1;
    end);

# Adjacency function for Latin square graphs.
BindGlobal("LatinSquareAdjacency",
    function(x, y)
        return x <> y and (x[1] = y[1] or x[2] = y[2]
                            or x[1]*x[2] = y[1]*y[2]);
    end);

# Adjacency function for adjacency lists.
BindGlobal("ListAdjacency",
    L -> function(x, y)
        return y in L[x];
    end);

# Point-line incidence
BindGlobal("PointLineIncidence", function(x, y)
    return x in y or y in x;
end);

# Adjacency for crooked graphy
BindGlobal("CrookedAdjacency",
    f -> function(x, y)
            return x <> y and x[3]+y[3] = f(x[1]+y[1])
                                    + (x[2]+y[2]+Z(2)^0)*(f(x[1]) + f(y[1]));
    end);
