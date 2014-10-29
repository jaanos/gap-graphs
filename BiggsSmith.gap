RequirePackage("grape");
I := [1,2,4,5,8,10];
G := Graph(Group(()), Cartesian([0..16], I),
    function(x,y) return x; end,
    function(x,y)
        if x = y then
            return false;
        fi;
        if x[1] = y[1] then
            return (x[2] in [5, 10] or y[2] in [5, 10]) and AbsoluteValue(x[2]-y[2]) in I;
        fi;
        return x[2] = y[2] and (not x[2] in [5, 10]) and (((x[1]-x[2]) mod 17) = y[1] or ((x[1]+x[2]) mod 17) = y[1]);
    end, true);
