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

# A bijective map from elements of F_q to integers from [0..q-1].
BindGlobal("FFEToInt", function(x, q)
    if IsZero(x) then
        return 0;
    else
        return LogFFE(x, Z(q))+1;
    fi;
end);

# A bijective map from integers from [0..q-1] to elements of F_q.
BindGlobal("IntToFFE", function(x, q)
    if x = 0 then
        return 0*Z(q);
    else
        return Z(q)^(x-1);
    fi;
end);

# The vector product of two vectors with 3 elements.
BindGlobal("VectorProduct", function(u, v)
    return [u[2]*v[3]-u[3]*v[2], u[3]*v[1]-u[1]*v[3], u[1]*v[2]-u[2]*v[1]];
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
BindGlobal("OnMatrices", function(q, d, e)
    local ij;
    ij := Cartesian([1..d], [1..e], [0..q-1]);
    return function(M, g)
        local ji, v, w, vw;
        ji := ij{OnTuples([1..d*e*q], g)};
        vw := List([1..d],
            i -> List([1..e],
                j -> ji[FFEToInt(M[i][j], q)+(j-1)*q+(i-1)*q*e+1]));
        v := List(vw, r -> List(r, x -> IntToFFE(x[3], q)));
        w := List(vw, r -> List(r, x -> x[2]));
        return List([1..d], i -> v[i]{List([1..e], j -> Position(w[i], j))});
    end;
end);

# Action of a wreath product on Hermitean matrices over a field.
BindGlobal("OnHermiteanMatrices", function(r, d)
    local c, ij;
    ij := Cartesian([1..d], [1..d], [0..r-1]);
    c := Conjugates(GF(r^2), GF(r), Z(r^2));
    return function(M, g)
        local ji, vw, F;
        ji := ij{OnTuples([1..d*d*r], g)};
        vw := List([1..d],
            i -> List([1..d],
                j -> ji[(j-1)*r+(i-1)*r*d+1]));
        F := function(i, j)
            if i = j then
                return IntToFFE(vw[i][j][3], r);
            elif i < j then
                return IntToFFE(vw[i][j][3], r) + c[1]*IntToFFE(vw[j][i][3], r);
            else
                return IntToFFE(vw[j][i][3], r) + c[2]*IntToFFE(vw[i][j][3], r);
            fi;
        end;
        return List([1..d],
            i -> List([1..d],
                j -> M[vw[i][j][1]][vw[i][j][2]] + F(i, j)));
    end;
end);

# Action of a matrix group on subspaces of a vector space over a finite field.
BindGlobal("OnSubspaces",
    V -> function(S, g)
        return Subspace(V, OnSubspacesByCanonicalBasis(Basis(S), g));
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

# Action on the vertices of additive symplectic covers of complete graphs.
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

# Action on the roots of E_8
BindGlobal("OnRoots", function(v, g)
    return (-1)^(10^g) * v{OnTuples([1..Length(v)], g)};
end);

# The field addition group as a permutation group.
BindGlobal("FieldAdditionPermutationGroup",
    q -> Group(List(Elements(Basis(GF(q))),
        g -> Permutation(g, GF(q), function(x, y) return x+y; end)))
);

# The field multiplication group as a permutation group.
BindGlobal("FieldMultiplicationPermutationGroup",
    q -> CyclicGroup(IsPermGroup, q-1)
);

# The field exponentiation group as a permutation group.
BindGlobal("FieldExponentiationPermutationGroup",
    q -> Action(Group(FrobeniusAutomorphism(GF(q))), GF(q))
);

# The group of even permutations of columns of a matrix
# as a permutation group over matrix elements.
BindGlobal("MatrixColumnEvenPermutationGroup", function(d, e)
    return Action(AlternatingGroup(e), Cartesian([1..d], [1..e]),
        function(x, g)
            return [x[1], OnPoints(x[2], g)];
        end);
end);

# The group of permutations of columns of a matrix
# as a permutation group over matrix elements.
BindGlobal("MatrixColumnPermutationGroup", function(d, e)
    return Action(SymmetricGroup(e), Cartesian([1..d], [1..e]),
        function(x, g)
            return [x[1], OnPoints(x[2], g)];
        end);
end);

# The group of simultaneous permutations of rows and columns
# of a square matrix as a permutation group over matrix elements.
BindGlobal("MatrixRowColumnPermutationGroup", function(d)
    return Action(SymmetricGroup(d), Cartesian([1..d], [1..d]), OnTuples);
end);

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

# Transforms a matrix over GF(r) to a Hermitean matrix over GF(r^2).
BindGlobal("ToHermitean", function(A, r)
    local c, n, Hermitize;
    c := Conjugates(GF(r^2), GF(r), Z(r^2));
    Hermitize := function(i, j)
        if i = j then
            return A[i][j];
        elif i < j then
            return A[i][j] + c[1]*A[j][i];
        else
            return A[j][i] + c[2]*A[i][j];
        fi;
    end;
    n := Size(A);
    return Immutable(List([1..n],
        i -> List([1..n], j -> Hermitize(i, j))));
end);

# Checks whether the sum of two subspaces are hyperbolic
# given a quadratic form.
BindGlobal("IsHyperbolic", function(Q)
    local s;
    s := (Size(BaseField(Q)) - 1)^2;
    return function(x, y)
        return Size(Filtered(x+y, z -> not IsZero(z*Q*z))) = s;
    end;
end);

# Checks whether two subspaces of F_{r^2} are orthogonal
# given a sesquilinear form.
BindGlobal("IsOrthogonal", function(Q, r)
    local F;
    F := x -> List(x, y -> y^r);
    return function(x, y)
        return ForAll(Cartesian(x, y), z -> IsZero(z[1]*Q*F(z[2])));
    end;
end);

# The subset of isotropic spaces with respect to the quadratic form Q
# of the collection V.
BindGlobal("IsotropicSpacesQuadraticForm",
    Q -> V -> Filtered(V, y -> ForAll(y, x -> IsZero(x*Q*x)))
);

# The subset of isotropic spaces with respect to the bilinear form Q
# of the collection V.
BindGlobal("IsotropicSpacesBilinearForm",
    Q -> V -> Filtered(V,
                y -> ForAll(Cartesian(y, y), x -> IsZero(x[1]*Q*x[2])))
);

# The subset of isotropic spaces with respect to the scalar product with
# conjugation map f of the collection V.
BindGlobal("IsotropicSpacesSesquilinearForm", function(Q, r)
    local O;
    O := IsOrthogonal(Q, r);
    return V -> Filtered(V, y -> O(y, y));
end);

# Adjacency function for Kneser-type graphs
BindGlobal("DisjointSets", function(x, y)
    return Intersection(x, y) = [];
end);

# Adjacency function for roots of E_8
BindGlobal("RootAdjacency", function(x, y)
    return x*y = 8;
end);

# Point-line incidence
BindGlobal("PointLineIncidence", function(x, y)
    return x in y or y in x;
end);

# Multiplication in Hall algebras
BindGlobal("HallMultiplication", function(p)
    local r;
    r := CoefficientsOfUnivariatePolynomial(p)[2];
    return function(x, y)
        if IsZero(y[2]) then
            return x*y[1];
        else
            return [x[1]*y[1] - x[2]/y[2]*Value(p, y[1]),
                    x[1]*y[2] - x[2]*(y[1] + r)];
        fi;
    end;
end);

# Multiplication in Dickson near-fields
BindGlobal("DicksonMultiplication",
    q -> function(x, y)
            if IsZero(y) then
                return 0*Z(q);
            else
                return x^(q^LogFFE(y, Z(q^2))) * y;
            fi;
        end);

# Right division in Dickson near-fields
BindGlobal("DicksonRightDivision",
    q -> function(x, y)
            if IsZero(x) then
                return 0*Z(q);
            else
                return (x / y)^(q^LogFFE(y, Z(q^2)));
            fi;
        end);

BindGlobal("ToExceptionalMatrix", function(q, F, B)
    return x -> F[IntVecFFE(Coefficients(B, x))*[q, 1] + 1];
end);

# Multiplication in exceptional near-fields
BindGlobal("ExceptionalMultiplication", function(q, F, B)
    local mat;
    mat := ToExceptionalMatrix(q, F, B);
    return function(x, y)
        local M;
        M := mat(x) * mat(y);
        return M[1]*B;
    end;
end);

# Right division in exceptional near-fields
BindGlobal("ExceptionalRightDivision", function(q, F, B)
    local mat;
    mat := ToExceptionalMatrix(q, F, B);
    return function(x, y)
        local M;
        M := mat(x) * mat(y)^-1;
        return M[1]*B;
    end;
end);

# Normalize a vector over a semifield given the semifield right division.
BindGlobal("NormalizeSemifieldVector",
    div -> function(v)
        local n;
        n := Filtered(v, x -> not IsZero(x))[1];
        return List(v, x -> div(x, n));
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
    local p1, p2, p3, act, F;
    p1 := Projection(dp, 1);
    p2 := Projection(dp, 2);
    p3 := Projection(dp, 3);
    act := OnSemifieldVectors(div);
    F := [Inverse, TransposedMat];
    return function(p, g)
        return [p[1]^Image(p2, g), act(List(p[2],
                        x -> x^(q^(2^Image(p3, g)))), F[p[1]](Image(p1, g)))];
    end;
end);
