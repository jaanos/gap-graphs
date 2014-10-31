RequirePackage("grape");

# The Doro graph with intersection array {12,10,3; 1,3,8} of nonisotropic
# 1-dimensional subspaces of F_4^4 with respect to a quadratic form and its
# associated bilinear form.
Doro4_Q := [[Z(4),1,1,1],[0,0,1,1],[0,0,0,1],[0,0,0,0]]*Z(4)^0;;
Doro4_V := Filtered(Elements(GF(4)^4), x -> x*Doro4_Q*x = Z(4)^0);;
Doro4 := Graph(Group(()), Doro4_V,
    function(x,y)
        return x;
    end,
    function(x,y)
        return x+y in Doro4_V;
    end, true);;

# The Doro graph with intersection array {10,6,4; 1,2,5} of nonisotropic
# 1-dimensional subspaces of F_5^4 with respect to a quadratic form and its
# associated bilinear form.
Doro5_Q := DiagonalMat([Z(5),1,1,1])*Z(5)^0;
Doro5_V := Filtered(List(Subspaces(GF(5)^4, 1), x -> Elements(x)[2]),
    z -> (z*Doro5_Q*z)^2 = Z(5)^0);;
Doro5 := Graph(Group(()), Doro5_V, function(x,y) return x; end,
    function(x,y)
        return x*Doro5_Q*y = 0*Z(5);
    end, true);;