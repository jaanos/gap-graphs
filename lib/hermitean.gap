RequirePackage("grape");

toHermitean := function(A, r)
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
end;

HFG := function(d, r)
    return Graph(Group(()), List(GF(r)^[d, d], x -> toHermitean(x, r)),
        function(x, y) return x; end,
        function(x, y)
            return RankMat(x-y) = 1;
        end, true);
end;