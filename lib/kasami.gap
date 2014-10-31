RequirePackage("grape");

Kasami := function(q, j, m)
    return Graph(Group(()), Elements(GF(q^(2*j+1))^2), function(x,y) return x; end,
        function(x,y)
            return x <> y and x[1]+y[1] = (x[2]+y[2])^(q^m+1);
        end, true);
end;
