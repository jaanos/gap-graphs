# De Caen, Mathon and Moorhouse's Preparata graph Pr(t, e)
InstallMethod(PreparataGraphCons,
    "as a spaces graph with full automorphism group", true,
    [IsSpacesGraph and FullAutomorphismGroup, IsInt, IsInt], 0,
    function(filter, t, e)
        local H;
        if t = 1 then
            H := HypercubeGraph(IsVectorGraph and FullAutomorphismGroup, 3);
            AssignVertexNames(H, List(H.names,
                                t -> ([t[1]+t[3], t[2]+t[3], t[3]-1])*Z(2)^0));
            return H;
        elif e mod t = 0 or t = 2 then
            TryNextMethod();
        else
            return PreparataGraphCons(IsSpacesGraph, t, e);
        fi;
    end);

InstallMethod(PreparataGraphCons, "as a spaces graph", true,
    [IsSpacesGraph, IsInt, IsInt], 0, function(filter, t, e)
        local q, s, F, K, dp;
        q := 2^(2*t-1);
        s := 2^e;
        F := GF(q);
        dp := DirectProduct(Group(Z(q)), FieldAdditionPermutationGroup(2),
                            FieldAdditionPermutationGroup(q),
                            FieldExponentiationPermutationGroup(q));
        return Graph(dp, Cartesian(F, GF(2), F), OnPreparata(q, s, dp),
                        CrookedAdjacency(GoldFunction(s)), true);
    end);

InstallMethod(PreparataGraphCons, "with full automorphism group", true,
    [FullAutomorphismGroup, IsInt, IsInt], 0, function(filter, t, e)
        return PreparataGraphCons(IsSpacesGraph
                                            and FullAutomorphismGroup, t, e);
    end);

InstallMethod(PreparataGraphCons, "default", true,
    [IsObject, IsInt, IsInt], 0, function(filter, t, e)
        return PreparataGraphCons(IsSpacesGraph, t, e);
    end);

BindGlobal("PreparataGraph", function(arg)
    local j, filt;
    if IsAFilter(arg[1]) then
        filt := arg[1];
        j := 2;
    else
        filt := IsObject;
        j := 1;
    fi;
    if Length(arg) = j then
        return PreparataGraphCons(filt, arg[j], 1);
    elif Length(arg) = j+1 then
        return PreparataGraphCons(filt, arg[j], arg[j+1]);
    elif Length(arg) = j+2 then
        return PreparataGraphCons(filt, arg[j], arg[j+1], arg[j+2]);
    else
        Error("usage: PreparataGraph( [<filter>, ]{<int>, <int>|<graph>}[, <int>] )");
    fi;
end);

# Quotient graph of the Preparata graph
BindGlobal("PreparataQuotientGraph", function(arg)
    local h, s, t, e, B, G, H, K, T, V;
    if Length(arg) < 2 then
        Error("at least two arguments expected");
        return fail;
    fi;
    h := arg[1];
    if IsGraph(arg[2]) then
        G := arg[2];
        if IsFFE(G.names[1][3]) then
            s := 1;
        else
            s := Size(G.names[1][3]);
        fi;
        t := Log2Int(2*s*G.order)/4;
    else
        t := arg[2];
        if Length(arg) > 2 then
            e := arg[3];
        else
            e := 1;
        fi;
        G := PreparataGraph(t, e);
    fi;
    B := BasisVectors(Basis(GF(2^(2*t-1))));
    if IsFFE(G.names[1][3]) then
        K := AdditiveGroup(B{[1..h]});
        V := List(G.names, x -> [x[1], x[2], x[3]+K]);
    else
        K := AdditiveGroup(B{[1..h+Log2Int(s)]});
        V := List(G.names, x -> [x[1], x[2], Elements(x[3])[1]+K]);
    fi;
    T := Set(List(V, x -> Positions(V, x)));
    H := Graph(Stabilizer(G.group, T, OnSetsSets), T, OnSets, function(x,y)
        return IsSubset(DistanceSet(G, 1, x), y);
    end, true);
    AssignVertexNames(H, List(T, x -> V[x[1]]));
    return H;
end);

# The coset graph of a Kasami code over an odd power extension
# of a binary field.
BindGlobal("KasamiGraph", function(i, j, m)
    local q, s, t, dp, G;
    q := 2^i;
    s := q^(2*j+1);
    t := q^m + 1;
    G := FieldAdditionPermutationGroup(s);
    dp := DirectProduct(G, G, FieldExponentiationPermutationGroup(s));
    return Graph(dp, Elements(GF(s)^2), OnKasami(s, s, dp),
        function(x, y)
            return x <> y and x[1]+y[1] = (x[2]+y[2])^t;
        end, true);
end);

# The coset graph of an extended Kasami code over an odd power extension
# of a binary field.
BindGlobal("ExtendedKasamiGraph", function(i, j, m)
    return ExtendedBipartiteDoubleGraph(KasamiGraph(i, j, m));
end);

# The coset graph of a Kasami code over a quadratic extension
# of a binary field.
BindGlobal("QuadraticKasamiGraph", function(i)
    local q, s, t, dp;
    q := 2^i;
    s := q^2;
    t := q + 1;
    dp := DirectProduct(FieldAdditionPermutationGroup(q),
            FieldAdditionPermutationGroup(s),
            FieldExponentiationPermutationGroup(s));
    return Graph(dp, Cartesian(GF(q), GF(s)), OnKasami(q, s, dp),
        function(x, y)
            return x <> y and x[1]+y[1] = (x[2]+y[2])^t;
        end, true);
end);

# The coset graph of an extended Kasami code over a quadratic extension
# of a binary field.
BindGlobal("ExtendedQuadraticKasamiGraph", function(i)
    return ExtendedBipartiteDoubleGraph(QuadraticKasamiGraph(i));
end);
