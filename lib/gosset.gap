RequirePackage("grape");

V := List(Combinations([1..8], 2),
    x -> List([1..8], function(i)
        if i in x then
            return 3;
        else
            return -1;
        fi;
    end));;
Append(V, -V);
G := Graph(Group(()), V, function(x, y) return x; end,
    function(x, y) return x*y = 8; end, true);;
