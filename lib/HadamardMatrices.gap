# The Sylvester Hadamard matrices of order 2^n
DeclareGlobalFunction("SylvesterHadamardMatrix");
InstallGlobalFunction(SylvesterHadamardMatrix, function(n)
    local E, H, H2, M, R, T, m, gens, gens2;
    H2 := [[1,1], [1,-1]];
    gens2 := [(1,2)(3,4)(6,8), (2,4)(5,6)(7,8)];
    if n < 0 then
        Error("nonnegative integer expected");
        return fail;
    elif n = 0 then
        return rec(matrix := ImmutableMatrix(Integers, [[1]]),
                   group := Group((1,2)(3,4)),
                   transpose := (1,3)(2,4));
    elif n = 1 then
        return rec(matrix := ImmutableMatrix(Integers, H2),
                   group := Group(gens2),
                   transpose := (1,5)(2,6)(3,7)(4,8));
    else
        m := Int(n/2);
        H := SylvesterHadamardMatrix(m);
        M := KroneckerProduct(H.matrix, H.matrix);
        T := PermList(Concatenation([2*2^n+1..4*2^n], [1..2*2^n]));
        if n mod 2 = 0 then
            R := Cartesian([1..2^m], [1..2^m]);
            return rec(matrix := ImmutableMatrix(Integers, M),
                   group := Action(Group(Union(List(GeneratorsOfGroup(H.group),
                                        g -> [DirectProductElement([g, ()]),
                                            DirectProductElement([(), g])]))),
                                Concatenation(Cartesian([1, -1], R),
                                                Cartesian([1, -1], R + 2^n)),
                                OnHadamardIndices([2^m,2^m])),
                   transpose := T);
        else
            gens := Union(Union(List(GeneratorsOfGroup(H.group),
                                g -> [DirectProductElement([g, (), ()]),
                                    DirectProductElement([(), g, ()])])),
                          List(gens2, g -> DirectProductElement([(), (), g])));
            R := Cartesian([1..2^m], [1..2^m], [1, 2]);
            return rec(matrix := ImmutableMatrix(Integers,
                                                 KroneckerProduct(M, H2)),
                   group := Action(Group(gens),
                                Concatenation(Cartesian([1, -1], R),
                                        Cartesian([1, -1], R + [2^n, 2^n, 4])),
                                OnHadamardIndices([2^m,2^m,2])),
                   transpose := T);
        fi;
    fi;
end);
