# Check whether we are given a filter (including IsObject)
BindGlobal("IsAFilter", x -> IsFilter(x) or x = IsObject);

# The canonical ordering of elements in F_q.
BindGlobal("FiniteFieldCanonicalOrdering",
    q -> Concatenation([0*Z(q)], List([0..q-2], i -> Z(q)^i)));

# A bijective map from elements of F_q to integers from [1..q].
BindGlobal("FFEToInt", function(x, q)
    if IsZero(x) then
        return 1;
    else
        return LogFFE(x, Z(q))+2;
    fi;
end);

# A bijective map from integers from [1..q] to elements of F_q.
BindGlobal("IntToFFE", function(x, q)
    if x = 1 then
        return 0*Z(q);
    else
        return Z(q)^(x-2);
    fi;
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

BindGlobal("ToExceptionalMatrix", function(q, F, B)
    return x -> F[IntVecFFE(Coefficients(B, x))*[q, 1] + 1];
end);

# Finds x and y in F_q such that x^2 + y^2 = z.
BindGlobal("AsSumOfSquares", function(z, q)
    local i, x, y;
    if IsZero(z) then
        return [z, z];
    elif q = 2 then
        return [0*Z(q), z];
    elif q mod 2 = 0 then
        return [0*Z(q), z^(1/2 mod (q-1))];
    elif IsOne(z^((q-1)/2)) then
        return [0*Z(q), Z(q)^(LogFFE(z, Z(q))/2)];
    fi;
    i := (q-1)/2;
    while i > 0 do
        x := Z(q)^i;
        y := z - x^2;
        if IsOne(y^((q-1)/2)) then
            return [x, Z(q)^(LogFFE(y, Z(q))/2)];
        fi;
        i := i-1;
    od;
    return fail;
end);

# Read a file into a list of lines.
BindGlobal("ReadLines", function(file)
    local f, ln, lns;
    f := InputTextFile(file);
    lns := [];
    ln := ReadLine(f);
    while ln <> fail do
        RemoveCharacters(ln, "\r\n");
        Add(lns, ln);
        ln := ReadLine(f);
    od;
    CloseStream(f);
    return lns;
end);
