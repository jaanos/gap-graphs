RequirePackage("grape");

Coxeter_P := List(Subspaces(GF(2)^3, 1), x -> List(x)[2]);
Coxeter_T := Filtered(Combinations(Coxeter_P, 3), x -> x[1] + x[2] <> x[3]);

# The Coxeter graph with intersection array {3,2,2,1; 1,1,1,2}
# as a graph of triangles of the Fano plane,
# with two triangles adjacent whenever they are disjoint.
Coxeter := Graph(Group(()), Coxeter_T, function(x,y) return x; end,
    function(x,y) return Length(Intersection(x, y)) = 0; end, true);;

Gewirtz_TT := Cartesian(Coxeter_T, GF(2));

Gewirtz_SignedAdjacency := function(x,y)
    return (x[2] = y[2] and Length(Intersection(x[1], y[1])) = 0)
        or (x[2] <> y[2] and x[1] = y[1])
        or (Length(Intersection(x[1], y[1])) = 1 and
            ((x[2] = 0*Z(2) and y[2] = Z(2)^0 and Sum(x[1]) in y[1])
            or (x[2] = Z(2)^0 and y[2] = 0*Z(2) and Sum(y[1]) in x[1])));
end;

# The Gewirtz graph with v = 56, k = 10, lm = 0, mu = 2
# constructed from two copies of the Coxeter graph.
Gewirtz := Graph(Group(()), Gewirtz_TT, function(x,y) return x; end,
    Gewirtz_SignedAdjacency, true);;
