RequirePackage("grape");

UnitaryNonisotropics := function(q)
    local P;
    P := Filtered(List(Subspaces(GF(q^2)^3, 1), z -> Elements(z)[2]),
        x -> Sum(List(x, y -> y^(q+1))) <> 0*Z(q));;
    return Graph(Group(()), P, function(x,y) return x; end,
        function(x,y)
            return Sum(List([1..3], i -> x[i]*y[i]^q)) = 0*Z(q);
    end, true);
end;