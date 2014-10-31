RequirePackage("grape");
# The Biggs-Smith graph with intersection array {3,2,2,2,1,1,1; 1,1,1,1,1,1,3}
# Constructed from 17 H shapes by connecting corresponding endpoints with
# cycles in four different manners.
BiggSmith_HShape := [1,2,4,5,8,10];
BiggsSmith := Graph(Group(()), Cartesian([0..16], BiggSmith_HShape),
    function(x,y) return x; end,
    function(x,y)
        if x = y then
            return false;
        fi;
        if x[1] = y[1] then
            return (x[2] in [5, 10] or y[2] in [5, 10]) and AbsoluteValue(x[2]-y[2]) in BiggSmith_HShape;
        fi;
        return x[2] = y[2] and (not x[2] in [5, 10]) and (((x[1]-x[2]) mod 17) = y[1] or ((x[1]+x[2]) mod 17) = y[1]);
    end, true);
