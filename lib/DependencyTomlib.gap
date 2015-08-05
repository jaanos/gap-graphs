if LoadPackage("tomlib") = true then
    # The Perkel graph with intersection array {6, 5, 2; 1, 1, 3}.
    InstallMethod(PerkelGraphCons, "as a conjugacy class graph", true,
        [IsConjugacyClassGraph], 1, function(filter)
            local G;
            G := PSL(2, 19);
            return Graph(G, Elements(First(ConjugacyClassesSubgroups(G),
                                    x -> Size(x) = 57 and Order(x[1]) = 60)),
                         ConjugateGroup, GroupIntersection(10), true);
        end);

    # The Biggs-Smith graph with intersection array
    # {3,2,2,2,1,1,1; 1,1,1,1,1,1,3}.
    InstallMethod(BiggsSmithGraphCons, "as a conjugacy class graph", true,
        [IsConjugacyClassGraph], 1, function(filter)
            local G;
            G := PSL(2, 17);
            return Graph(G, Elements(First(ConjugacyClassesSubgroups(G),
                                    x -> Size(x) = 102 and Order(x[1]) = 24)),
                         ConjugateGroup, GroupIntersection(8), true);
        end);
fi;
