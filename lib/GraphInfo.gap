# Checks whether a graph is an antipodal cover.
BindGlobal("IsAntipodalCover", function(G)
    local d, k, ia, i, ci, cj;
    if not IsSimpleGraph(G) then
        Error("not a simple graph");
        return fail;
    fi;
    if not IsConnectedGraph(G) then
        Error("not a connected graph");
        return fail;
    fi;
    d := Diameter(G);
    ia := GlobalParameters(G);
    k := ia[1][3];
    if k = -1 or ia[d+1][1] <> k then
        return false;
    fi;
    for i in [1..Length(G.representatives)] do
        ci := DistanceSet(G, 1, DistanceSet(G, [0, d], G.representatives[i]));
        cj := Union(List(Adjacency(G, G.representatives[i]),
                x -> DistanceSet(G, [0, d], x)));
        if ci <> cj then
            return false;
        fi;
    od;
    return true;
end);

# Covering index of an antipodal cover.
BindGlobal("AntipodalCoveringIndex", function(G)
    if not IsAntipodalCover(G) then
        Error("not an antipodal cover");
        return fail;
    fi;
    return Length(DistanceSet(G, [0, Diameter(G)], 1));
end);

# Parameters of a generalized polygon.
# Does not check whether G actually is a generalized polygon
# (for d = 2 this is not guaranteed).
BindGlobal("GeneralizedPolygonParameters", function(G)
    local d, s, t, ia;
    if not IsSimpleGraph(G) then
        Error("not a simple graph");
        return fail;
    fi;
    if not IsConnectedGraph(G) then
        Error("not a connected graph");
        return fail;
    fi;
    d := Diameter(G);
    ia := GlobalParameters(G);
    t := ia[d+1][1] - 1;
    s := ia[1][3] / ia[d+1][1];
    if not IsInt(s) or ForAny(ia{[2..d]}, x -> x[1] <> 1 or x[3] <> s*t) then
        return fail;
    fi;
    return [2*d, s, t];
end);

# Check whether two vertices of a generalized quadrangle with parameter t
# are either equal, or make a regular pair.
BindGlobal("IsRegularPair", function(G, x, y, t)
    return x = y or y in Adjacency(G, x) or
        Size(Intersection(List(Intersection(Adjacency(G, x), Adjacency(G, y)),
            z -> Adjacency(G, z)))) = t+1;
end);

# Find regular points in a generalized quadrangle.
BindGlobal("RegularPoints", function(G)
    local l, t, P, orb, par;
    par := GeneralizedPolygonParameters(G);
    if par = fail then
        return fail;
    fi;
    t := par[3];
    l := [1..Length(G.representatives)];
    P := [1..G.order];
    orb := List(l, i -> OrbitsDomain(Stabilizer(G.group,
                                                    G.representatives[i]), P));
    return G.representatives{Filtered(l, i -> ForAll(orb[i][1],
                        x -> IsRegularPair(G, G.representatives[i], x, t)))};
end);
