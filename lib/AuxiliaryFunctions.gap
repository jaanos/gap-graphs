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
