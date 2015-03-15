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

# The incidence graph of an ordinary Hughes plane.
# It is distance-regular when q is an odd prime power.
BindGlobal("HughesPlaneIncidenceGraph", function(q)
    local c, A, B, P, dp, th, mul, rdiv;
    c := CoefficientsOfUnivariatePolynomial(DefiningPolynomial(GF(GF(q), 3)));
    A := [[-c[3], Z(q)^0, 0*Z(q)],
          [-c[2], 0*Z(q), Z(q)^0],
          [-c[1], 0*Z(q), 0*Z(q)]];
    mul := DicksonMultiplication(q);
    rdiv := DicksonRightDivision(q);
    P := Filtered(GF(q^2)^3,
            x -> not IsZero(x) and IsOne(Filtered(x, y -> not IsZero(y))[1]));
    th := Z(q^2)^((q+1)/2);
    B := Basis(GF(q^2), [Z(q)^0, th]);
    dp := DirectProduct(Group(A), SymmetricGroup(2), SymmetricGroup(2));
    return Graph(dp, Cartesian([1, 2], P),
            OnHughesPlane(q, rdiv, dp),
            function(x, y)
                local z;
                z := TransposedMat(List(y[2], w -> Coefficients(B, w)));
                return x[1] <> y[1] and IsZero(x[2]*z[1] + mul(th, x[2]*z[2]));
            end, true);;
end);

# The incidence graph of an exceptional Hughes plane. The first argument q must
# be 5, 7, 11, 23, 29, or 59. As there are two exceptional near-fields of order
# 11^2, a second parameter equal to 1 or 2 may be supplied when q = 11.
BindGlobal("ExceptionalHughesPlaneIncidenceGraph", function(arg)
    local c, n, q, A, B, F, G, P, dp, th, mul, rdiv, gens;
    if Length(arg) < 1 then
        Error("at least one argument expected");
        return fail;
    fi;
    q := arg[1];
    if Length(arg) > 1 then
        n := arg[2];
    else
        n := 1;
    fi;
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
    c := CoefficientsOfUnivariatePolynomial(DefiningPolynomial(GF(GF(q), 3)));
    A := [[-c[3], Z(q)^0, 0*Z(q)],
          [-c[2], 0*Z(q), Z(q)^0],
          [-c[1], 0*Z(q), 0*Z(q)]];
    P := Filtered(GF(q^2)^3,
            x -> not IsZero(x) and IsOne(Filtered(x, y -> not IsZero(y))[1]));
    th := Z(q^2)^((q+1)/2);
    B := Basis(GF(q^2), [Z(q)^0, th]);
    mul := ExceptionalMultiplication(q, F, B);
    rdiv := ExceptionalRightDivision(q, F, B);
    dp := DirectProduct(Group(A), SymmetricGroup(2), SymmetricGroup(1));
    return Graph(dp, Cartesian([1, 2], P),
            OnHughesPlane(q, rdiv, dp),
            function(x, y)
                local z;
                z := TransposedMat(List(y[2], w -> Coefficients(B, w)));
                return x[1] <> y[1] and IsZero(x[2]*z[1] + mul(th, x[2]*z[2]));
            end, true);;
end);
