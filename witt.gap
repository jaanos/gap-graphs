RequirePackage("grape");

B := [
[1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,1,1,1,1],
[0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,1,0,0,1,1],
[0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,0,0,1],
[0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,1,1,0,1],
[0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,1,1,1],
[0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,1,1,1],
[0,0,0,0,0,0,1,0,0,0,0,0,1,1,1,0,0,1,1,0,1,1,0,0],
[0,0,0,0,0,0,0,1,0,0,0,0,1,1,1,1,0,0,0,1,0,1,1,0],
[0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,1,0,1,0,1,0,1,0],
[0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,1,1,1,1,0,1,0,0],
[0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,1,1,0,1,1,0,1,0],
[0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,1]]*Z(2)^0;;

V := VectorSpace(GF(2), B);
P := Filtered(V, x -> WeightVecFFE(x) = 8);;
Q := Filtered(P, x -> x[1] = 0*Z(2));;
R := Filtered(Q, x -> x[2] = 0*Z(2));;

W24 := Graph(Group(()), P, function(x,y) return x; end,
    function(x,y)
        return WeightVecFFE(x+y) = 16;
    end, true);;
W23 := Graph(Group(()), Q, function(x,y) return x; end,
    function(x,y)
        return WeightVecFFE(x+y) = 16;
    end, true);;
W22 := Graph(Group(()), R, function(x,y) return x; end,
    function(x,y)
        return WeightVecFFE(x+y) = 16;
    end, true);;
