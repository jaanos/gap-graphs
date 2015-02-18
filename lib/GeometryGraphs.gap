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

# The incidence graph of a Hughes plane.
# It is distance-regular when q is an odd prime power.
BindGlobal("HughesPlaneIncidenceGraph", function(q)
    local c, A, G, L, P, V, act, dp, mul, rdiv;
    c := CoefficientsOfUnivariatePolynomial(DefiningPolynomial(GF(GF(q), 3)));
    A := [[-c[3], Z(q)^0, 0*Z(q)],
          [-c[2], 0*Z(q), Z(q)^0],
          [-c[1], 0*Z(q), 0*Z(q)]];
    V := GF(q^2)^3;
    mul := DicksonMultiplication(q);
    rdiv := DicksonRightDivision(q);
    P := Filtered(V,
            x -> not IsZero(x) and IsOne(Filtered(x, y -> not IsZero(y))[1]));
    L := List(Filtered(GF(q^2), w -> IsOne(w) or not w in GF(q)),
            t -> Set(Filtered(P, x -> IsZero(x[1] + mul(t, x[2]) + x[3]))));
    G := Group(A);
    act := OnSetsSemifieldVectors(rdiv);
    dp := DirectProduct(G, Group((1, 2)));
    return Graph(dp, Concatenation(P, Union(List(L,
                                        l -> Set(List(G, g -> act(l, g)))))),
                OnHughesPlane(q, rdiv, dp), PointLineIncidence, true);
end);
