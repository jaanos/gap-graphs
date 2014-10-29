RequirePackage("grape");

P := List(Subspaces(GF(2)^3, 1), x -> List(x)[2]);
T := Filtered(Combinations(P, 3), x -> x[1] + x[2] <> x[3]);

Coxeter := Graph(Group(()), T, function(x,y) return x; end,
    function(x,y) return Length(Intersection(x, y)) = 0; end, true);;

TT := Cartesian(T, GF(2));

SignedAdjacency := function(x,y)
    return (x[2] = y[2] and Length(Intersection(x[1], y[1])) = 0)
        or (x[2] <> y[2] and x[1] = y[1])
        or (Length(Intersection(x[1], y[1])) = 1 and
            ((x[2] = 0*Z(2) and y[2] = Z(2)^0 and Sum(x[1]) in y[1])
            or (x[2] = Z(2)^0 and y[2] = 0*Z(2) and Sum(y[1]) in x[1])));
end;

Gewirtz := Graph(Group(()), TT, function(x,y) return x; end,
    SignedAdjacency, true);;


O := [[1,8,11,13,15],[1,9,10,12,14],[2,4,7,14,15],[2,5,6,12,13],[3,4,6,10,11],[3,5,7,8,9]];
S := [[[1,2,3], [4,8,12], [5,11,14], [7,10,13], [6,9,15]],
      [[1,2,3], [5,10,15], [4,9,13], [6,8,14], [7,11,12]],
      [[1,4,5], [2,8,10], [3,13,14], [7,11,12], [6,9,15]],
      [[1,4,5], [2,9,11], [3,12,15], [6,8,14], [7,10,13]],
      [[1,6,7], [2,8,10], [3,12,15], [4,9,13], [5,11,14]],
      [[1,6,7], [2,9,11], [3,13,14], [4,8,12], [5,10,15]]];
p := [[[1,2,3],[2,3,5]], [[1,2,4],[2,4,6]], [[1,2,5],[1,3,6]], [[1,2,6],[1,4,5]],
      [[1,3,4],[1,5,6]], [[1,3,5],[1,2,4]], [[1,3,6],[3,4,6]], [[1,4,5],[3,4,5]],
      [[1,4,6],[1,2,3]], [[1,5,6],[2,5,6]], [[2,3,4],[1,3,4]], [[2,3,5],[4,5,6]],
      [[2,3,6],[1,2,6]], [[2,4,5],[1,2,5]], [[2,4,6],[3,5,6]], [[2,5,6],[2,3,4]],
      [[3,4,5],[2,3,6]], [[3,4,6],[2,4,5]], [[3,5,6],[1,3,5]], [[4,5,6],[1,4,6]]];
      
G2 := Graph(Group(()), Union(Cartesian([0], O, S), List(p, x -> Concatenation([1], x))),
    function(x,y) return x; end, function(x,y)
        local o, s;
        if x[1] = y[1] then
            if x[1] = 0 then
                o := Intersection(x[2], y[2]);
                s := Intersection(x[3], y[3]);
                return Length(o) = 1 and Length(s) = 1 and o[1] in s[1];
            else
                return Length(Intersection(x[2], y[2])) = 0;
            fi;
        else
            if x[1] = 0 then
                s := x;
                x := y;
                y := s;
            fi;
            return y[2] in O{x[2]} and y[3] in S{x[3]};
        fi;
    end, true);;
    
xy := [1, 2];
uv := [7, 8];
Nxy := DistanceSet(Gewirtz, 1, xy);
Nuv := DistanceSet(Gewirtz, 1, uv);
M := Intersection(Nxy, Nuv);
Mxy := Difference(Nxy, M);
Muv := Difference(Nuv, M);
wz := Filtered(Mxy, x -> Length(Intersection(Mxy, Adjacency(Gewirtz, x))) = 1);
st := Filtered(Muv, x -> Length(Intersection(Muv, Adjacency(Gewirtz, x))) = 1);
Mx := Intersection(Adjacency(Gewirtz, xy[1]), Union(M, wz));
My := Intersection(Adjacency(Gewirtz, xy[2]), Union(M, wz));
Dx := Difference(Union(List(Mx, x -> Union(List(Filtered(Mx, z -> z <> x),
        y -> Intersection(Adjacency(Gewirtz, x), Adjacency(Gewirtz, y)))))),
        Union(xy, uv));
Dy := Difference(Union(List(My, x -> Union(List(Filtered(My, z -> z <> x),
        y -> Intersection(Adjacency(Gewirtz, x), Adjacency(Gewirtz, y)))))),
        Union(xy, uv));
S := Difference(Union([xy, uv, Mxy, Muv, Dx, Dy]), Union(wz, st));

Syl2 := InducedSubgraph(Gewirtz, S);

l := [[0*Z(2), 0*Z(2), Z(2)^0], [0*Z(2), Z(2)^0, 0*Z(2)], [0*Z(2), Z(2)^0, Z(2)^0]];
TS := Filtered(TT, x -> Length(Intersection(x[1], l)) in [2, 2-IntFFE(x[2])]);

Sylvester := Graph(Group(()), TS, function(x,y) return x; end,
    SignedAdjacency, true);;
