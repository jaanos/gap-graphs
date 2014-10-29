#!/usr/bin/gap -o 2G

RequirePackage("grape");
V := GF(2)^8;
es := Elements(Subspaces(V, 4));;
M := [
    [ 0, 1, 1, 1, 1, 1, 1, 1 ],
    [ 1, 0, 1, 1, 1, 1, 1, 1 ],
    [ 1, 1, 0, 1, 1, 1, 1, 1 ],
    [ 1, 1, 1, 0, 1, 1, 1, 1 ],
    [ 1, 1, 1, 0, 0, 1, 1, 1 ],
    [ 1, 1, 0, 1, 1, 0, 1, 1 ],
    [ 1, 0, 1, 1, 1, 1, 0, 1 ],
    [ 0, 1, 1, 1, 1, 1, 1, 0 ] ];
iso := Filtered(es, e -> not Z(2)^0 in List(Elements(e), x -> x*M*x));;
Length(iso);
# 270

G := Graph(Group(()), iso, function(x,y) return x; end,
    function(x,y) return Dimension(Intersection(x,y)) = 3; end, true);;
GlobalParameters(G);
# [ [ 0, 0, 15 ], [ 1, 0, 14 ], [ 3, 0, 12 ], [ 7, 0, 8 ], [ 15, 0, 0 ] ]

SG := Graph(Group(()), DistanceSet(G, 4, 1), function(x,y) return x; end,
    function(x,y) return Distance(G, x, y) = 2; end, true);;
CG := ComplementGraph(SG);;

LG := EdgeGraph(SG);; # geometrija z d=4, g=3
LG := EdgeGraph(CG);; # geometrija z d=4, g=3
