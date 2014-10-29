RequirePackage("grape");

Q := [[Z(4),1,1,1],[0,0,1,1],[0,0,0,1],[0,0,0,0]]*Z(4)^0;;
V := Filtered(Elements(GF(4)^4), x -> x*Q*x = Z(4)^0);;
G := Graph(Group(()), V,
    function(x,y)
        return x;
    end,
    function(x,y)
        return x+y in V;
    end, true);;
GlobalParameters(G);