RequirePackage("grape");

OnGraphs := function(A, g)
    return OnTuplesSets(A, g){OnTuples([1..Size(A)], g^-1)};
end;
    
OnTuplesGraphs := function(T, g)
    return List(T, G -> OnGraphs(G, g));
end;

R := List([(1,2,4,5,6,7), (1,4,5,6,7,3), (1,2,5,4,6,7), (1,4,2,5,6,7),
           (1,3,5,4,6,7), (1,4,3,5,6,2), (1,3,5,4,6,2)],
    r -> List(List([1..7],
        x -> Set([x^r, x^(r^4)])),
            function(s)
                if Length(s) = 1 then
                    return [];
                else
                    return s;
                fi;
            end));

S7 := SymmetricGroup(7);
C := Orbits(S7, [[[2,5],[3,6],[1,4],[2,5],[3,6],[1,4],[]]], OnGraphs)[1];;
L := Union(Filtered(Orbits(S7, List(C, c -> [C[1], c]), OnTuplesGraphs),
                    o -> o[1][2] in R));;

G := Graph(Group(()), C,
    function(x,y)
        return x;
    end,
    function(x,y)
        return [x,y] in L;
    end, true);;
GlobalParameters(G);
