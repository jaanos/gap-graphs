# The incidence graph of a Desarguesian projective plane.
BindGlobal("DesarguesianPlaneIncidenceGraph", function(q)
    local V, dp;
    V := GF(q)^3;
    dp := DirectProduct(GL(3, q), SymmetricGroup(2));
    return Graph(dp, Union(List([1,2], d -> Subspaces(V, d))),
                 OnProjectivePlane(V, dp), function(x, y)
                    return x <> y and Intersection(x, y) in [x, y];
                 end, true);
end);

# The incidence graph of a Hall plane.
BindGlobal("HallPlaneIncidenceGraph", function(q)
    local c, p, H, L, P, dp, mul;
    H := GF(q)^2;
    p := DefiningPolynomial(GF(GF(q), 2));
    c := CoefficientsOfUnivariatePolynomial(p);
    mul := HallMultiplication(p);
    P := Union([[]], List(H, x -> [x]), Cartesian(H, H));
    L := Union([Union([[]], List(H, z -> [z]))],
                List(H, x -> Union([[]], List(H, z -> [x, z]))),
                List(Cartesian(H, H), w -> Union([[w[1]]],
                    List(H, z -> [z, mul(z, w[1]) + w[2]]))));
    dp := DirectProduct(Concatenation(ListWithIdenticalEntries(5,
                                        FieldAdditionPermutationGroup(q)),
                        [FieldMultiplicationPermutationGroup(q),
                         Group([[c[2], Z(q)^0], [-c[1], 0*Z(q)]])]));
    return Graph(dp, Union(P, L), OnHallPlane(q, dp),
                PointLineIncidence, true);
end);

# The incidence graph of a Hughes plane. If the second parameter n is zero or
# unspecified, a Dickson semifield of order q^2 is used and the first parameter
# q must be an odd prime power. Otherwise, an exceptional near-field of order
# q^2 is used, so q must be 5, 7, 11, 23, 29, or 59. As there are two
# exceptional near-fields of order 11^2, setting n to 1 or 2 chooses one of
# these semifields when q = 11.
BindGlobal("HughesPlaneIncidenceGraph", function(arg)
    local c, n, q, A, B, F, G, P, dp, th, mul, rdiv, gens;
    if Length(arg) < 1 then
        Error("at least one argument expected");
        return fail;
    fi;
    q := arg[1];
    if Length(arg) > 1 then
        n := arg[2];
    else
        n := 0;
    fi;
    th := Z(q^2)^((q+1)/2);
    B := Basis(GF(q^2), [Z(q)^0, th]);
    if n = 0 then
        mul := DicksonMultiplication(q);
        rdiv := DicksonRightDivision(q);
    else
        gens := [[[0, -1], [1, 0]]];
        if q = 5 then
            gens[2] := [[1, -2], [-1, -2]];
        elif q = 7 then
            gens[2] := [[1, 3], [-1, -2]];
        elif q = 11 and n = 1 then
            gens[2] := [[1, 5], [-5, -2]];
            gens[3] := [[4, 0], [0, 4]];
        elif q = 11 and n = 2 then
            gens[2] := [[2, 4], [1, -3]];
        elif q = 23 then
            gens[2] := [[1, -6], [12, -2]];
            gens[3] := [[2, 0], [0, 2]];
        elif q = 29 then
            gens[2] := [[1, -7], [-12, -2]];
            gens[3] := [[16, 0], [0, 16]];
        elif q = 59 then
            gens[2] := [[9, 15], [-10, -10]];
            gens[3] := [[4, 0], [0, 4]];
        else
            Error("the specified exceptional near-field does not exist");
            return fail;
        fi;
        G := Group(gens * Z(q)^0);
        F := Union([NullMat(2, 2, GF(q))], Elements(G));
        SortBy(F, x -> IntVecFFE(x[1]));
        mul := ExceptionalMultiplication(q, F, B);
        rdiv := ExceptionalRightDivision(q, F, B);
    fi;
    c := CoefficientsOfUnivariatePolynomial(DefiningPolynomial(GF(GF(q), 3)));
    A := [[-c[3], Z(q)^0, 0*Z(q)],
          [-c[2], 0*Z(q), Z(q)^0],
          [-c[1], 0*Z(q), 0*Z(q)]];
    P := Filtered(GF(q^2)^3,
            x -> not IsZero(x) and IsOne(Filtered(x, y -> not IsZero(y))[1]));
    dp := DirectProduct(Group(A), SymmetricGroup(2));
    return Graph(dp, Cartesian([1, 2], P),
            OnHughesPlane(q, rdiv, dp),
            function(x, y)
                local z;
                z := TransposedMat(List(y[2], w -> Coefficients(B, w)));
                return x[1] <> y[1] and IsZero(x[2]*z[1] + mul(th, x[2]*z[2]));
            end, true);;
end);

# The collinearity graph of the generalized quadrangle Q(d, q)
# of order (q, q^{d-3}).
BindGlobal("GeneralizedQuadrangleQ", function(d, q)
    return PolarGraphO(4-d, d+1, q);
end);

# The collinearity graph of the generalized quadrangle H(d, r^2)
# of order (r^2, r^{d-5/2}).
BindGlobal("GeneralizedQuadrangleH", function(d, r)
    return PolarGraphU(d+1, r);
end);

# The collinearity graph of the generalized quadrangle W(q) of order (q, q).
BindGlobal("GeneralizedQuadrangleW", q -> PolarGraphSp(4, q));

# The collinearity graph of the generalized quadrangle T_d(O) of order
# (q, q^{d-1}) derived from the projective space PG(d+1, q) containing the
# oval or ovoid O in a hyperplane.
BindGlobal("GeneralizedQuadrangleT", function(arg)
    local d, o, q, G, H, L, O, P, V, dp;
    if Length(arg) < 3 then
        Error("at least three arguments expected");
        return fail;
    fi;
    d := arg[1];
    q := arg[2];
    o := arg[3];
    if Length(arg) > 3 then
        G := arg[4];
    else
        G := Group(());
    fi;
    V := GF(q)^(d+2);
    H := Subspace(V, BasisVectors(CanonicalBasis(V)){[2..d+2]}, "basis");
    O := List(o, x -> Subspace(V,
            [Concatenation([0*Z(q)], BasisVectors(Basis(x))[1])], "basis"));
    P := Union(Filtered(Subspaces(V, 1), x -> not IsSubset(H, x)),
                Filtered(Subspaces(V, d+1),
                    x -> Length(Filtered(O, y -> IsSubset(x, y))) = 1),
                [Subspace(V, [], "basis")]);
    return Graph(Group(()), P, function(x, y) return x; end,
        function(x, y)
            local dx, dy, xy, yx;
            dx := Dimension(x);
            dy := Dimension(y);
            xy := x+y;
            yx := Intersection(x, y);
            return (dx = 0 and dy = d+1) or (dx = d+1 and dy = 0) or
                (dx = 1 and dy = 1 and ForAny(O, z -> IsSubset(xy, z))) or
                (dx = 1 and dy = d+1 and ForAny(O, z -> IsSubset(y, x+z))) or
                (dx = d+1 and dy = 1 and ForAny(O, z -> IsSubset(x, y+z))) or
                (dx = d+1 and dy = d+1 and x <> y and
                                            ForAny(O, z -> IsSubset(yx, z)));
        end, true);
end);

# The collinearity graph of the generalized quadrangle P(G, z) of
# order (s-1, s+1) derived by removing the neighbourhood of a regular point z
# of a generalized quadrangle G of order (s, s).
BindGlobal("GeneralizedQuadrangleP", function(G, z)
    local H;
    H := Graph(Stabilizer(G.group, z), DistanceSet(G, 2, z), OnPoints,
            function(x, y)
                local c;
                c := Intersection(Adjacency(G, x), Adjacency(G, y));
                return not z in c and (y in Adjacency(G, x)
                                    or ForAll(c, w -> z in Adjacency(G, w)));
            end, true);
    return H;
end);

# The collinearity graph of the generalized quadrangle AS(q)
# of order (q-1, q+1).
BindGlobal("GeneralizedQuadrangleAS",
    q -> GeneralizedQuadrangleP(GeneralizedQuadrangleW(q), 1));

# The incidence graph of a projective plane read from a file as on
#   http://www.uwyo.edu/moorhouse/pub/planes/       or
#   http://www.uwyo.edu/moorhouse/pub/genpoly/
# The optional second parameter is a file containing generators of the
# automorphism group.
BindGlobal("IncidenceGraphFromFile", function(arg)
    local l, m, n, G, L, fst, lst, lns;
    if Length(arg) < 1 then
        Error("at least one argument expected");
        return fail;
    fi;
    lns := ReadLines(arg[1]);
    n := Length(lns);
    L := List(lns, l -> List(SplitString(l, " "), x -> Int(x)+n+1));
    if Length(arg) > 1 then
        lns := ReadLines(arg[2]);
        m := Length(lns);
        fst := Filtered([1..m], i -> IntChar(lns[i][1]) <> 32);
        l := Length(fst);
        lst := List(fst{[2..l]}, i -> i-1);
        Add(lst, m);
        G := Group(List([1..l], i -> PermList(List(SplitString(
                    JoinStringsWithSeparator(lns{[fst[i]..lst[i]]}, ""),
                    " "), x -> Int(x)+1))));
    else
        G := Group(());
    fi;
    return Graph(G, [1..2*n], OnPoints, function(x, y)
            return (x <= n and y in L[x]) or (y <= n and x in L[y]);
        end, true);
end);

# The collinearity graph of a projective plane read from a file as on
#   http://www.uwyo.edu/moorhouse/pub/genpoly/
# The optional second parameter is a file containing generators of the
# automorphism group.
BindGlobal("CollinearityGraphFromFile", function(arg)
    local l, m, n, G, L, fst, lst, lns;
    if Length(arg) < 1 then
        Error("at least one argument expected");
        return fail;
    fi;
    lns := ReadLines(arg[1]);
    n := Length(lns);
    L := List(lns, l -> List(SplitString(l, " "), x -> Int(x)));
    if Length(arg) > 1 then
        lns := ReadLines(arg[2]);
        m := Length(lns);
        fst := Filtered([1..m], i -> IntChar(lns[i][1]) <> 32);
        l := Length(fst);
        lst := List(fst{[2..l]}, i -> i-1);
        Add(lst, m);
        G := Group(List([1..l], i -> PermList(Filtered(List(SplitString(
                    JoinStringsWithSeparator(lns{[fst[i]..lst[i]]}, ""),
                    " "), x -> Int(x)+1), y -> y <= n))));
    else
        G := Group(());
    fi;
    return Graph(G, [1..n], OnPoints, function(x, y)
            return Length(Intersection(L[x], L[y])) = 1;
        end, true);
end);

