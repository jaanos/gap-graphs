# Action on a finite field element.
BindGlobal("OnFFE", q -> function(x, g)
    return IntToFFE(FFEToInt(x, q)^g, q);
end);

# Action of a group on a signed point.
BindGlobal("OnSignedPoints", function(dp, signs)
    local p1, p2;
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    return function(e, g)
        return [OnPoints(e[1], Image(p1, g)),
                Permuted(signs, Image(p2, g))[Position(signs, e[2])]];
    end;
end);

# Action of a product group on vertices of a product graph.
BindGlobal("OnProduct", function(n, dp)
    local pr;
    pr := List([1..n], i -> Projection(dp, i));
    return function(e, g)
        return List([1..n],
            i -> OnPoints(e[i], Image(pr[i], g)));
    end;
end);

# Action of a product group on vertices of a sum graph.
BindGlobal("OnSum", dp -> function(e, g)
        return [e[1], OnPoints(e[2], Image(Projection(dp, e[1]), g))];
    end
);

# Action of a product group on the Cayley table of its first factor.
BindGlobal("OnLatinSquare", function(dp)
    local p1, p2, p3, p4, p5;
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    p3 := Projection(dp, 3);
    p4 := Projection(dp, 4);
    p5 := Projection(dp, 5);
    return function(e, g)
        local g3, g5, s, t;
        g3 := Image(p3, g);
        g5 := Image(p5, g);
        s := SignPerm(g5);
        t := Permuted([Image(p1, g) * e[1]^g3, e[2]^g3 * Image(p2, g)],
                        Image(p4, g));
        t[3] := (t[1]*t[2])^-1;
        return List(Permuted(t, g5){[1,2]}, x -> x^s);
    end;
end);

# Action on the vertices of the Johnson graph J(2d, d).
BindGlobal("OnJohnson", function(n, dp)
    local p1, p2, F;
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    F := [s -> s, s -> Difference([1..n], s)];
    return function(s, g)
        return F[1^Image(p2, g)](OnSets(s, Image(p1, g)));
    end;
end);

# Action on the vertices of the Chang graphs.
BindGlobal("OnChang", function(dp)
    local p1, p2;
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    return function(s, g)
        local g2, t;
        g2 := Image(p2, g);
        s := OnSets(s, Image(p1, g) * g2);
        if Order(g2) = 2 then
            t := First([[1..4], [5..8]], x -> IsSubset(x, s));
            if t <> fail then
                s := Difference(t, s);
            fi;
        fi;
        return s;
    end;
end);

# Action of a wreath product on vectors over a ring.
BindGlobal("OnWreathProduct", function(d, e, wp)
    local p, q;
    p := Projection(wp);
    q := [0..d-1]*e;
    return function(v, g)
        return Permuted(List([0..d-1], i -> (v[i+1]+i*e)^g), Image(p, g)) - q;
    end;
end);

# Action on pairs of vectors with 3 elements.
BindGlobal("OnVectorPairs", function(q, dp, wp)
    local A, p, p1, p2, p3;
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    p3 := Projection(dp, 3);
    p := Projection(wp);
    A := OnWreathProduct(6, q, wp);
    return function(M, g)
        local u, v, w, g1, g2;
        g1 := Image(p1, g);
        g2 := Image(p2, g);
        u := List(Concatenation(M), x -> FFEToInt(x, q)^Image(p3, g));
        v := Permuted(List(u, x -> IntToFFE(x, q)), Image(p, g1));
        w := List(A(u, g1), x -> IntToFFE(x, q));
        return [Determinant(g2)^-1 * g2 *
                    (w{[1..3]} + VectorProduct(v{[4..6]}, w{[4..6]})),
                TransposedMat(g2)^-1 * w{[4..6]}];
    end;
end);

# Action on matrices over a field.
BindGlobal("OnMatrices", function(q, d, e, dp)
    local F, p1, p2, pm;
    F := OnFFE(q);
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    pm := List([0..d-1], i -> List([1..e], j -> Projection(dp, e*i+j+2)));
    return function(M, g)
        return Image(p1, g)*M*Image(p2, g) +
            List(pm, r -> List(r, p -> F(0*Z(q), Image(p, g))));
    end;
end);

# Action on Hermitean matrices over a field.
BindGlobal("OnHermiteanMatrices", function(r, d, dp)
    local C, F, K, p1, pm;
    F := OnFFE(r);
    K := Elements(GF(r^2));
    C := List(K, x -> Conjugates(GF(r^2), GF(r), x)[2]);
    p1 := Projection(dp, 1);
    pm := List([0..d-1], i -> List([1..d], j -> Projection(dp, d*i+j+1)));
    return function(M, g)
        local H;
        H := Image(p1, g);
        return List(H, v -> List(v, x -> C[Position(K, x)]))*M*TransposedMat(H)
            + ToHermitean(List(pm, v -> List(v, p -> F(0*Z(r), Image(p, g)))),
                            r);
    end;
end);

# Action of a matrix group on subspaces of a vector space over a finite field.
BindGlobal("OnSubspaces",
    V -> function(S, g)
        return Subspace(V, OnSubspacesByCanonicalBasis(Basis(S), g));
    end);

# Action of a matrix group on a set of subspaces of a vector space over a
# finite field.
BindGlobal("OnSetsSubspaces", function(V)
    local A;
    A := OnSubspaces(V);
    return function(L, g)
        return Set(List(L, S -> A(S, g)));
    end;
end);

# Action on the vertices of doubled Odd graphs.
BindGlobal("OnDoubledOdd", function(n, dp)
    local F, p1, p2;
    F := [x -> x, x -> Difference([1..n], x)];
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    return function(s, g)
        return F[1^Image(p2, g)](OnSets(s, Image(p1, g)));
    end;
end);

# Action on the vertices of doubled Grassmann graphs.
BindGlobal("OnDoubledGrassmann", function(V, dp)
    local A, F, p1, p2;
    F := [x -> x, OrthogonalSpaceInFullRowSpace];
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    A := OnSubspaces(V);
    return function(S, g)
        return F[1^Image(p2, g)](A(S, Image(p1, g)));
    end;
end);

# Action on the vertices of Paley graphs.
BindGlobal("OnPaley", function(q, dp)
    local F, p1, p2, p3;
    F := OnFFE(q);
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    p3 := Projection(dp, 3);
    return function(x, g)
        return (F(x, Image(p3, g) * Image(p1, g))) * Z(q)^((q-1)^Image(p2, g));
    end;
end);

# Action on the vertices of Preparata graphs.
BindGlobal("OnPreparata", function(q, s, dp)
    local F, p1, p2, p3, p4;
    F := OnFFE(q);
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    p3 := Projection(dp, 3);
    p4 := Projection(dp, 4);
    return function(t, g)
        local z;
        z := Z(q)^((q-1)^Image(p1, g));
        return List([z*t[1], t[2] + 2^Image(p2, g),
                     F(t[3]*z^(s+1), Image(p3, g))], x -> F(x, Image(p4, g)));
    end;
end);

# Action on the vertices of Kasami graphs.
BindGlobal("OnKasami", function(q, s, dp)
    local Fq, Fs, p1, p2, p3;
    Fq := OnFFE(q);
    Fs := OnFFE(s);
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    p3 := Projection(dp, 3);
    return function(v, g)
        return List([v[1] + Fq(0*Z(q), Image(p1, g)), Fs(v[2], Image(p2, g))],
                    x -> Fs(x, Image(p3, g)));
    end;
end);

# Action on vertices of additive symplectic covers of complete graphs.
BindGlobal("OnAdditiveSymplecticCover", function(q, r, m, B, K, dp)
    local F, p1, p2, p3, p4, p5, pr;
    F := OnFFE(q);
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    p3 := Projection(dp, 3);
    p4 := Projection(dp, 4);
    p5 := Projection(dp, 5);
    pr := List([1..m], i -> Projection(dp, i+5));
    return function(p, g)
        local h, k, y, z, g4;
        g4 := Image(p4, g);
        z := List([1..m], i -> F(0*Z(q), Image(pr[i], g)));
        h := Z(q)^((q-1)^Image(p3, g));
        y := List(Concatenation(h*p[2]{[1..m/2]}, p[2]{[m/2+1..m]}),
                    x -> F(x, g4));
        k := Z(r)^((r-1)^Image(p2, g));
        return [K + k^2 * (F(h*Elements(p[1])[1], g4*Image(p5, g)) + y*B*z),
                k * (y*Image(p1, g) + z)];
    end;
end);

# Action on vertices of multiplicative symplectic covers of complete graphs.
BindGlobal("OnMultiplicativeSymplecticCover", function(q, dp)
    local F, p1, p2, p3;
    F := OnFFE(q);
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    p3 := Projection(dp, 3);
    return function(s, g)
        local g2, g3;
        g2 := Image(p2, g);
        g3 := Image(p3, g);
        return OnSets(Set(List(s, p -> List(Permuted([1,2], g2),
                                            i -> F(p[i], g3)))),
                      Image(p1, g));
    end;
end);

# Action on the vertices of affine polar graphs.
BindGlobal("OnAffine", function(q, d, dp)
    local F, p1, pr;
    F := OnFFE(q);
    p1 := Projection(dp, 1);
    pr := List([1..d], i -> Projection(dp, i+1));
    return function(v, g)
        return v^Image(p1, g) + List([1..d], i -> F(0*Z(q), Image(pr[i], g)));;
    end;
end);

# Action on the roots of E_8
BindGlobal("OnRoots", function(v, g)
    return (-1)^(10^g) * v{OnTuples([1..Length(v)], g)};
end);

# Action of a matrix group on normalized vectors over a semifield.
BindGlobal("OnSemifieldVectors", function(div)
    local norm;
    norm := NormalizeSemifieldVector(div);
    return function(s, g)
        return norm(s*g);
    end;
end);

# Action of a matrix group on sets of normalized vectors over a semifield.
BindGlobal("OnSetsSemifieldVectors", function(div)
    local norm;
    norm := NormalizeSemifieldVector(div);
    return function(S, g)
        return Set(List(S*g, norm));
    end;
end);

# Action on points and lines of Desarguesian projective planes.
BindGlobal("OnProjectivePlane", function(V, dp)
    local F, P, p1, p2;
    F := OnSubspaces(V);
    P := [x -> x, OrthogonalSpaceInFullRowSpace];
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    return function(S, g)
        return P[2^Image(p2, g)](F(S, Image(p1, g)));
    end;
end);

# Action on points or lines of a point-line geometry.
BindGlobal("OnPointsOrLines", function(act, line)
    return function(x, g)
        if line(x) then
            return Set(List(x, p -> act(p, g)));
        else
            return act(x, g);
        fi;
    end;
end);

# Action on points and lines of Hall planes.
BindGlobal("OnHallPlane", function(q, dp)
    local p5, p6, p7, p8, pr, A, F;
    pr := List([0,1], i -> List([1,2], j -> Projection(dp, 2*i+j)));
    p5 := Projection(dp, 5);
    p6 := Projection(dp, 6);
    p7 := Projection(dp, 7);
    F := OnFFE(q);
    A := function(p, g)
        local s, M;
        M := Image(p7, g);
        p := List(p, z -> [z[1] + F(0*Z(q), Image(p5, g))*z[2],
                            z[2]*(Z(q)^((q-1)^Image(p6, g)))]);
        if p = [] then
            if IsZero(M[2][1]) then
                return [];
            else
                return [[M[2][2]/M[1][2], 0*Z(q)]];
            fi;
        elif Length(p) = 1 then
            if not IsZero(p[1][2]) then
                return p;
            fi;
            s := M[1][1] + p[1][1]*M[1][2];
            if IsZero(s) then
                return [];
            else
                return [[(M[2][1] + p[1][1]*M[2][2])/s, 0*Z(q)]];
            fi;
        else
            return M*p + List(pr, r -> List(r, e -> F(0*Z(q), Image(e, g))));
        fi;
    end;
    return OnPointsOrLines(A, x -> Length(x) > 2);
end);

# Action on points and lines of Hughes planes.
BindGlobal("OnHughesPlane", function(q, div, dp)
    local p1, p2, act, F;
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    act := OnSemifieldVectors(div);
    F := [Inverse, TransposedMat];
    return function(p, g)
        return [p[1]^Image(p2, g), act(p[2], F[p[1]](Image(p1, g)))];
    end;
end);
