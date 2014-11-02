# Action of a product group on vertices of a product graph.
BindGlobal("OnProduct", function(n, dp)
    return function(e, g)
        return List([1..n],
            i -> OnPoints(e[i], Image(Projection(dp, i), g)));
    end;
end);

# Action of a product group on vertices of a sum graph.
BindGlobal("OnSum", dp -> function(e, g)
        return [e[1], OnPoints(e[2], Image(Projection(dp, e[1]), g))];
    end
);

# Action of a group on its multiplication table.
BindGlobal("OnLatinSquare", function(e, g)
    return [e[1] * g, g^-1 * e[2]];
end);

# Action of a wreath product on vectors over a ring.
BindGlobal("OnZmodnZVectors", function(d, e)
    local ij;
    ij := Cartesian([1..d], [0..e-1]);
    return function(vec, g)
        local ji, v, w, vw;
        ji := ij{OnTuples([1..d*e], g)};
        vw := List([1..d], i -> ji[Int(vec[i])+(i-1)*e+1]);
        v := List(vw, x -> ZmodnZObj(x[2], e));
        w := List(vw, x -> x[1]);
        return v{List([1..d], i -> Position(w, i))};
    end;
end);

# Transforms a matrix over GF(r) to a Hermitean matrix over GF(r^2).
BindGlobal("ToHermitean", function(A, r)
    local H, c, x, y, i, j, n;
    H := MutableCopyMat(A);
    n := Size(H);
    c := Conjugates(GF(r^2), GF(r), Z(r^2));
    for i in [1..n] do
        for j in [i+1..n] do
            x := H[i][j];
            y := H[j][i];
            H[i][j] := x + c[1]*y;
            H[j][i] := x + c[2]*y;
        od;
    od;
    return H;
end);
