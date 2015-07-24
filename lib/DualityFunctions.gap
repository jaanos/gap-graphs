# Default function for dual points.
BindGlobal("DefaultDualityFunction", x -> x);

# Default function for primal points.
BindGlobal("DefaultPrimalityFunction", x -> Intersection(x)[1]);

# Duality function for bipartite doubles.
BindGlobal("BipartiteDoubleDualityFunction",
    f -> x -> [x[1][1], f(List(x, y -> y[2]))]
);

# Duality function for Grassmann graphs.
BindGlobal("GrassmannDualityFunction", function(x)
    local y;
    y := Intersection(x);
    if Dimension(y) = 0 then
        return Sum(x);
    else
        return y;
    fi;
end);

# Check whether function for dual and primal points exist,
# and add them if they do not.
BindGlobal("CheckDualityFunctions", function(G)
    if not "duality" in RecNames(G) then
        G.duality := DefaultDualityFunction;
    fi;
    if not "primality" in RecNames(G) then
        G.primality := DefaultPrimalityFunction;
    fi;
end);
