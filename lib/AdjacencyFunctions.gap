# Adjacency function for Kneser-type graphs.
BindGlobal("DisjointSets", function(x, y)
    return Intersection(x, y) = [];
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

# Adjacency function for adjacency lists.
BindGlobal("ListAdjacency",
    L -> function(x, y)
        return y in L[x];
    end);

# Point-line incidence
BindGlobal("PointLineIncidence", function(x, y)
    return x in y or y in x;
end);
