RequirePackage("grape");

# bilinear forms graph: vozlišča so matrike,
# dve vozlišči sta sosednji, ko je njuna razlika ranga 1
GB := AdjFunGraph(Elements(GF(2)^[2, 3]),
    function(x,y) return RankMat(x-y) = 1; end);;

# OA(8, 3): vozlišča so vnosi v latinskem kvadratu (Z_8)^2,
# dve vozlišči sta sosednji, če sta v isti vrstici, v istem stolpcu,
# ali pa imata isto vrednost
GZ8 := AdjFunGraph(Elements(ZmodnZ(8)^2), function(x,y) return x <> y and
        (x[1] = y[1] or x[2] = y[2] or x[1]+x[2] = y[1]+y[2]); end);;

# OA(8, 3) za latinski kvadrat (Z_2 × Z_4)^2
A := Elements(AbelianGroup([2, 4]));;
A2 := Cartesian(A, A);;
GZ24 := AdjFunGraph(A2, function(x,y) return x <> y and
        (x[1] = y[1] or x[2] = y[2] or x[1]*x[2] = y[1]*y[2]); end);;

# OA(8, 3) za latinski kvadrat (Dih(8))^2
D := Elements(DihedralGroup(8));;
D2 := Cartesian(D, D);;
GD8 := AdjFunGraph(D2, function(x,y) return x <> y and
        (x[1] = y[1] or x[2] = y[2] or x[1]*x[2] = y[1]*y[2]); end);;

# OA(8, 3) za latinski kvadrat (Q_8)^2
Q := Elements(SmallGroup(8, 4));;
Q2 := Cartesian(Q, Q);;
GQ8 := AdjFunGraph(Q2, function(x,y) return x <> y and
        (x[1] = y[1] or x[2] = y[2] or x[1]*x[2] = y[1]*y[2]); end);;

# iz [21, 6; 8:21, 12:42]_2 kode
M := [
    [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    [1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    [1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0],
    [1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0],
    [1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0],
    [1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1]
]*Z(2);

# vozlišča so kodne besede,
# dve vozlišči sta sosednji, ko je razdalja med njima minimalna
GC := AdjFunGraph(Elements(Subspace(GF(2)^21, M)),
    function(x, y) return WeightVecFFE(x+y) = 8; end);;
lGC := DistanceSetInduced(GC, 1, 1);;
