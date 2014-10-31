RequirePackage("grape");
S := PSL(2,19);
V := Elements(ConjugacyClassesSubgroups(S)[16]);;
G := Graph(Group(()), V, function(x,y) return x; end,
    function(x,y)
        return x <> y and Order(Intersection(x,y)) mod 5 = 0;
    end, true);;
GlobalParameters(G);