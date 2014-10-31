RequirePackage("grape");

M := [[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]];
V := GF(4)^4;
P := Elements(Subspaces(V, 1));
L := Filtered(Elements(Subspaces(V, 2)), S -> Length(Filtered(Elements(S), x -> x*M*Elements(S)[2] = 0*Z(2))) = 16);
G := Graph(Group(()), P, function(x,y) return x; end, function(x,y) return x <> y and Length(Filtered(L, l -> IsSubset(l, x) and IsSubset(l, y))) > 0; end, true);;
GlobalParameters(G);

co := [];
K55 := [];
for i in [1..85] do
    for j in [i+1..85] do
        if Distance(G, i, j) = 1 then
            continue;
        fi;
        mu := Intersection(Adjacency(G, i), Adjacency(G, j));
        if mu in co then
            continue;
        fi;
        um := Intersection(List(mu, x -> Adjacency(G, x)));
        Add(co, um);
        Add(co, mu);
        Add(K55, [um, mu]);
    od;
od;

GG := Graph(Group(()), K55, function(x,y) return x; end, function(x,y)
    return (Length(Intersection(x[1], y[1])) = 1 and Length(Intersection(x[2], y[2])) = 1)
        or (Length(Intersection(x[1], y[2])) = 1 and Length(Intersection(x[2], y[1])) = 1);
    end, true);;
GlobalParameters(GG);
