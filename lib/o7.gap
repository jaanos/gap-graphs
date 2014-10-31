RequirePackage("grape");

V := GF(3)^7;
P := Filtered(List(Elements(Subspaces(V, 1)), v -> Elements(v)[2]), u -> u*u = Z(3));;

G := Graph(Group(()), P, function(x,y) return x; end, function(x,y) return x*y = 0*Z(3); end, true);;
GlobalParameters(G);
