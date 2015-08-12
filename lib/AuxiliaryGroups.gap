# The field addition group as a permutation group.
BindGlobal("FieldAdditionPermutationGroup",
    q -> Group(List(Elements(Basis(GF(q))),
        g -> Permutation(g, FiniteFieldCanonicalOrdering(q), \+)))
);

# The field exponentiation group as a permutation group.
BindGlobal("FieldExponentiationPermutationGroup",
    q -> Action(Group(FrobeniusAutomorphism(GF(q))),
                FiniteFieldCanonicalOrdering(q))
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

# The wreath product of two symmetric groups.
BindGlobal("WreathProductSymmetricGroups", function(m, n)
    local g;
    if not (IsInt(m) and IsInt(n) and m >= 0 and n >= 0) then
        Error("m or n not a nonnegative integer");
        return fail;
    fi;
    if m <= 1 then
        return SymmetricGroup(n);
    elif n <= 1 then
        return SymmetricGroup(m);
    fi;
    g := Concatenation(GeneratorsOfGroup(SymmetricGroup(m)),
                        [PermList(Concatenation([m+1..2*m], [1..m]))]);
    if n > 2 then
        Add(g, PermList(Concatenation([m+1..n*m], [1..m])));
    fi;
    return Group(g);
end);

# The automorphism group of a (hyper)conic
# in a Desarguesian projective plane of order q.
BindGlobal("HyperconicAutomorphismGroup", q -> Group([
    [[-Z(q)^0, 0*Z(q), 0*Z(q)],
     [ 0*Z(q), 0*Z(q), Z(q)^0],
     [ 0*Z(q), Z(q)^0, 0*Z(q)]],
    [[-Z(q)^0, 0*Z(q), -2*Z(q)^0],
     [ Z(q)^0, Z(q)^0,   Z(q)^0],
     [ 0*Z(q), 0*Z(q),   Z(q)^0]],
    [[-Z(q),   0*Z(q), -2*Z(q)^0],
     [ Z(q),   Z(q)^2,   Z(q)^0],
     [ 0*Z(q), 0*Z(q),   Z(q)^0]],
    Z(q)*IdentityMat(3, GF(q))
]));

BindGlobal("OrthogonalNegationGroup", function(arg)
    local H, M, S, d, e, h, q, s, t;
    if Length(arg) < 2 then
        Error("at least two arguments expected");
        return fail;
    elif Length(arg) = 2 then
        e := 0;
        d := arg[1];
        q := arg[2];
    else
        e := arg[1];
        d := arg[2];
        q := arg[3];
    fi;
    H := [[Z(q)^0, 0*Z(q)],
          [0*Z(q), -Z(q)^0]];
    if e = 0 or q mod 2 = 0 or (d = 2 and e = -1) then
        M := IdentityMat(d, Z(q));
    elif d = 2 then
        M := H;
    elif q mod 4 = 1 then
        M := Z(q)^((q-1)/4) * IdentityMat(d, Z(q));
        M[1][1] := Z(q)^0;
        M[2][2] := -Z(q)^0;
    else
        t := AsSumOfSquares(-Z(q)^0, q);
        h := d/2;
        S := [[t[1], t[2]],
              [t[2], -t[1]]];
        if (e = 1) = (d mod 4 = 0) then
            s := Z(q)^((q+1)/4);
            M := BlockMatrix(Concatenation([[1, 1, H],
                                            [2, 2, [[0*Z(q), s^-1],
                                                    [s, 0*Z(q)]]]],
                          List([3..h], i -> [i, i, S])), h, h);
        else
            t := AsSumOfSquares(-Z(q)^0, q);
            M := BlockMatrix(Concatenation([[1, 1, H]],
                          List([2..h], i -> [i, i, S])), h, h);
        fi;
    fi;
    return Group(M);
end);