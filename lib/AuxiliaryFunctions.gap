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
    if q mod 2 = 0 then
        return [Z(q)^0, (z-Z(q)^0)^(1/2 mod (q-1))];
    fi;
    i := (q-1)/2;
    while i > 0 do
        x := Z(q)^(2*i);
        y := z - x;
        if IsOne(y^((q-1)/2)) then
            return [x, y];
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
