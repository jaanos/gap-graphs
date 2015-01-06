# De Caen, Mathon and Moorhouse's Preparata graph Pr(t, e)
BindGlobal("PreparataGraph", function(arg)
    local t, e, q, s, F, K, dp;
    if Length(arg) < 1 then
        Error("at least one argument expected");
        return fail;
    fi;
    t := arg[1];
    if Length(arg) > 1 then
        e := arg[2];
    else
        e := 1;
    fi;
    q := 2^(2*t-1);
    s := 2^e;
    F := GF(q);
    dp := DirectProduct(FieldMultiplicationPermutationGroup(q),
                        FieldAdditionPermutationGroup(2),
                        FieldAdditionPermutationGroup(q),
                        FieldExponentiationPermutationGroup(q));
    return Graph(dp, Cartesian(F, GF(2), F), OnPreparata(q, s, dp),
        function(x,y)
            return x <> y
                and x[3]+y[3] = x[1]^s * y[1] + x[1] * y[1]^s + (x[2]+y[2])*(x[1]^(s+1) + y[1]^(s+1));
        end, true);
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
