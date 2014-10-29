RequirePackage("grape");

FlagGraph := pg -> Graph(Group(()),
    Union(List([1..Length(pg)], x -> List(pg[x], y -> [x, y]))),
    function(x,y) return x; end,
    function(x,y) return x <> y and (x[1] = y[1] or x[2] = y[2]); end, true);

IncGraph := pg -> Graph(Group(()),
    Union(List([1..91], x -> [0, x]), List(Union(pg), x -> [1, x])),
    function(x,y) return x; end,
    function(x,y)
        local p, l;
        if x[1] = y[1] then
            return false;
        fi;
        if x[1] = 0 then
            l := x[2];
            p := y[2];
        else
            p := x[2];
            l := y[2];
        fi;
        return p in pg[l];
    end, true);
