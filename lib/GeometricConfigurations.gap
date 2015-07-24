# A hyperconic in a Desarguesian projective plane of order 2^s.
BindGlobal("HyperconicProjectivePlane",
    s -> rec(points := Set(List(Union([[0,1,0], [1,0,0]]*Z(2^s)^0,
                                    List(GF(2^s), t -> [t, t^2, 1]*Z(2^s)^0)),
                                x -> Subspace(GF(2^s)^3, [x], "basis"))),
            group := HyperconicAutomorphismGroup(2^s), d := 2, q := 2^s)
);

# A classical conic in a Desarguesian projective plane of order q.
BindGlobal("ClassicalConicProjectivePlane",
    q -> rec(points := Set(List(Union([[0,1,0]]*Z(q)^0, List(GF(q),
            t -> [t, t^2, 1]*Z(q)^0)), x -> Subspace(GF(q)^3, [x], "basis"))),
            group := HyperconicAutomorphismGroup(q), d := 2, q := q)
);

# A nonclassical conic in a Desarguesian projective plane of order q.
BindGlobal("NonclassicalConicProjectivePlane",
    q -> rec(points := Set(List(Union([[1,0,0]]*Z(q)^0, List(GF(q),
            t -> [t, t^2, 1]*Z(q)^0)), x -> Subspace(GF(q)^3, [x], "basis"))),
            group := Group([DiagonalMat([Z(q), Z(q)^2, Z(q)^0]),
                            Z(q)*IdentityMat(3, GF(q))]),
            d := 2, q := q)
);

# An elliptic quadric in a d-dimensional projective space of order q.
BindGlobal("EllipticQuadricProjectiveSpace", function(q)
    local G, V;
    G := GO(-1, 4, q);
    V := GF(q)^4;
    return rec(points := IsotropicSpacesQuadraticForm(
                    InvariantQuadraticForm(G).matrix)(Subspaces(V, 1)),
                group := G, d := 3, q := q);
end);
