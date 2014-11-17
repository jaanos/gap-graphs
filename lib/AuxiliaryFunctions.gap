# Action of a group on a signed point.
BindGlobal("OnSignedPoints", function(dp, signs)
    return function(e, g)
        return [OnPoints(e[1], Image(Projection(dp, 1), g)),
            Permuted(signs, Image(Projection(dp, 2), g))[Position(signs, e[2])]];
    end;
end);

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

# Action of a product group on the multiplication table of its factors.
BindGlobal("OnLatinSquare", dp -> function(e, g)
        return [Image(Projection(dp, 1), g) * e[1],
            e[2] * Image(Projection(dp, 2), g)];
    end
);

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

BindGlobal("OnSubspaces",
    V -> function(S, g)
        return Subspace(V, OnSubspacesByCanonicalBasis(Basis(S), g));
    end
);

# The field addition group as a permutation group.
BindGlobal("FieldAdditionPermutationGroup",
    q -> Group(List(Elements(Basis(GF(q))),
        g -> Permutation(g, GF(q), function(x, y) return x+y; end)))
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

# The subset of isotropic spaces with respect to the quadratic form Q
# of the collection V.
BindGlobal("IsotropicSpacesQuadraticForm",
    Q -> V -> Filtered(V, y -> Size(Filtered(Elements(y),
            x -> not IsZero(x*Q*x))) = 0)
);

# The subset of isotropic spaces with respect to the bilinear form Q
# of the collection V.
BindGlobal("IsotropicSpacesBilinearForm",
    Q -> V -> Filtered(V,
                y -> Size(Filtered(Cartesian(Elements(y), Elements(y)),
                    x -> not IsZero(x[1]*Q*x[2]))) = 0)
);

# The subset of isotropic spaces with respect to the scalar product with
# conjugation map f of the collection V.
BindGlobal("IsotropicSpacesSesquilinearForm", function(Q, r)
    local F;
    F := x -> List(x, y -> y^r);
    return V -> Filtered(V,
                y -> Size(Filtered(Cartesian(Elements(y), Elements(y)),
                    x -> not IsZero(x[1]*Q*F(x[2])))) = 0);
end);
