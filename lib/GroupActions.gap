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

# Action of a product group on the multiplication table of its factors.
BindGlobal("OnLatinSquare", function(dp)
    local p1, p2;
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    return function(e, g)
        return [Image(p1, g) * e[1], e[2] * Image(p2, g)];
    end;
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

# Action of a wreath product on pairs of vectors with 3 elements.
BindGlobal("OnVectorPairs", function(q)
    local ij;
    ij := Cartesian([1..2], [1..3], [0..q-1]);
    return function(M, g)
        local ji, v, w, vw;
        ji := ij{OnTuples([1..6*q], g)};
        vw := List([1..2],
            i -> List([1..3],
                j -> ji[FFEToInt(M[i][j], q)+(j-1)*q+(i-1)*q*3+1]));
        v := List(vw, r -> List(r, x -> IntToFFE(x[3], q)));
        w := List(vw, r -> List(r, x -> x[2]));
        v := List([1..2], i -> v[i]{List([1..3], j -> Position(w[i], j))});
        w := List([1..2], i -> M[i]{List([1..3], j -> Position(w[i], j))});
        return [v[1] + VectorProduct(w[2], v[2]), v[2]];
    end;
end);

# Action of a wreath product on matrices over a field.
BindGlobal("OnMatrices", function(q, d, e, dp)
    local F, N, p1, p2, pm;
    F := Elements(GF(q));
    N := Position(F, 0*Z(q));
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    pm := List([0..d-1], i -> List([1..e], j -> Projection(dp, e*i+j+2)));
    return function(M, g)
        return Image(p1, g)*M*Image(p2, g) +
            List(pm, r -> List(r, p -> F[N^Image(p, g)]));
    end;
end);

# Action of a wreath product on Hermitean matrices over a field.
BindGlobal("OnHermiteanMatrices", function(r, d, dp)
    local C, F, K, N, p1, pm;
    F := Elements(GF(r));
    N := Position(F, 0*Z(r));
    K := Elements(GF(r^2));
    C := List(K, x -> Conjugates(GF(r^2), GF(r), x)[2]);
    p1 := Projection(dp, 1);
    pm := List([0..d-1], i -> List([1..d], j -> Projection(dp, d*i+j+1)));
    return function(M, g)
        local H;
        H := Image(p1, g);
        return List(H, r -> List(r, x -> C[Position(K, x)]))*M*TransposedMat(H)
            + ToHermitean(List(pm, r -> List(r, p -> F[N^Image(p, g)])), r);
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
    local F, N, p1, p2;
    F := Elements(GF(q));
    N := Position(F, 0*Z(q));
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    return function(x, g)
        return (x + F[N^Image(p1, g)]) * Z(q)^((q-1)^Image(p2, g));
    end;
end);

# Action on the vertices of Preparata graphs.
BindGlobal("OnPreparata", function(q, s, dp)
    local F, N, p1, p2, p3, p4;
    F := Elements(GF(q));
    N := Position(F, 0*Z(q));
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    p3 := Projection(dp, 3);
    p4 := Projection(dp, 4);
    return function(t, g)
        local w, z;
        z := Z(q)^((q-1)^Image(p1, g));
        w := [z*t[1], t[2] + 2^Image(p2, g),
              t[3]*z^(s+1) + F[N^Image(p3, g)]];
        return List(w, x -> F[Position(F, x)^Image(p4, g)]);
    end;
end);

# Action on the vertices of Kasami graphs.
BindGlobal("OnKasami", function(q, s, dp)
    local Fq, Fs, Nq, Ns, p1, p2, p3;
    Fq := Elements(GF(q));
    Nq := Position(Fq, 0*Z(q));
    Fs := Elements(GF(s));
    Ns := Position(Fs, 0*Z(s));
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    p3 := Projection(dp, 3);
    return function(v, g)
        local w;
        w := [v[1] + Fq[Nq^Image(p1, g)],
              v[2] + Fs[Ns^Image(p2, g)]];
        return List(w, x -> Fs[Position(Fs, x)^Image(p3, g)]);
    end;
end);

# Action on vertices of additive symplectic covers of complete graphs.
BindGlobal("OnAdditiveSymplecticCover", function(q, m, B, dp)
    local F, N, p1, p2, pr;
    F := Elements(GF(q));
    N := Position(F, 0*Z(q));
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    pr := List([1..m], i -> Projection(dp, i+2));
    return function(p, g)
        local z;
        z := List([1..m], i -> F[N^Image(pr[i], g)]);
        return [p[1] + F[N^Image(p2, g)] + p[2]*B*z,
                p[2]*Image(p1, g) + z];
    end;
end);

# Action on vertices of multiplicative symplectic covers of complete graphs.
BindGlobal("OnMultiplicativeSymplecticCover", function(q, dp)
    local F, p1, p2, p3;
    F := Elements(GF(q));
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    p3 := Projection(dp, 3);
    return function(s, g)
        local g2, g3;
        g2 := Image(p2, g);
        g3 := Image(p3, g);
        return OnSets(Set(List(s, p -> List(Permuted([1,2], g2),
                                            i -> F[Position(F, p[i])^g3]))),
                      Image(p1, g));
    end;
end);

# Action on the vertices of affine polar graphs.
BindGlobal("OnAffine", function(q, d, dp)
    local F, N, p1, pr;
    F := Elements(GF(q));
    N := Position(F, 0*Z(q));
    p1 := Projection(dp, 1);
    pr := List([1..d], i -> Projection(dp, i+1));
    return function(v, g)
        return v^Image(p1, g) + List([1..d], i -> F[N^Image(pr[i], g)]);;
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
    local p5, p6, p7, p8, pr, A, F, N;
    pr := List([0,1], i -> List([1,2], j -> Projection(dp, 2*i+j)));
    p5 := Projection(dp, 5);
    p6 := Projection(dp, 6);
    p7 := Projection(dp, 7);
    F := Elements(GF(q));
    N := Position(F, 0*Z(q));
    A := function(p, g)
        local s, M;
        M := Image(p7, g);
        p := List(p, z -> [z[1] + F[N^Image(p5, g)]*z[2],
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
            return M*p + List(pr, r -> List(r, e -> F[N^Image(e, g)]));
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
