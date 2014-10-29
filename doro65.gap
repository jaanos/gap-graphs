RequirePackage("grape");

Q := DiagonalMat([Z(5),1,1,1])*Z(5)^0;
V := Filtered(List(Subspaces(GF(5)^4, 1), x -> Elements(x)[2]),
    z -> (z*Q*z)^2 = Z(5)^0);;
G := Graph(Group(()), V, function(x,y) return x; end,
    function(x,y)
        return x*Q*y = 0*Z(5);
    end, true);;
GlobalParameters(G);