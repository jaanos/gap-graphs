if LoadPackage("loops") then
    # Action of a product group on the Cayley table of a loop.
    BindGlobal("OnLatinSquareFromLoop", function(dp)
        local p1, p2, p3, p4, p5;
        p1 := Projection(dp, 1);
        p2 := Projection(dp, 2);
        p3 := Projection(dp, 3);
        p4 := Projection(dp, 4);
        p5 := Projection(dp, 5);
        return function(e, g)
            local g2, g4;
            g2 := Image(p2, g);
            g4 := Image(p4, g);
            return [Image(p1, g) * e[1]^g4 * g2,
                    g2^-1 * e[2]^g4 * Image(p3, g)];
        end;
    end);
        
    # Latin square graphs from quasigroup.
    InstallMethod(LatinSquareGraphCons, "for loops", true,
    [IsObject, IsQuasigroup, IsBool], 0, function(filter, Q, invt)
            local CL, CQ, G, L, e, g1, g2;
            e := Elements(Q);
            CQ := CayleyTable(Q);
            CL := NormalizedQuasigroupTable(CQ);
            L := LoopByCayleyTable(CL);
            G := LatinSquareGraphCons(filter, L, invt);
            g1 := PermList(CQ[1])^-1;
            g2 := PermList(CQ{[1..Length(CQ)]}[1^g1])^-1;
            AssignVertexNames(G, List(G.names, t -> [e[Position(L, t[1])^g2],
                                                    e[Position(L, t[2])^g1]]));
            return G;
        end);
    
    # Latin square graphs from loops.
    InstallMethod(LatinSquareGraphCons, "for loops", true,
    [IsObject, IsLoop, IsBool], 0, function(filter, L, invt)
            local A, G, dp, vcs;
            if IsAbelian(L) then
                A := Group((1,2));
            else
                A := Group(());
            fi;
            G := l -> GroupByGenerators(GeneratorsOfLoop(l), One(l));
            dp := DirectProduct(G(LeftNucleus(L)), G(MiddleNucleus(L)),
                                G(RightNucleus(L)), AutomorphismGroup(L), A);
            if invt then
                vcs := Cartesian(L, L);
            else
                vcs := [[One(L), One(L)]];
            fi;
            return Graph(dp, vcs, OnLatinSquareFromLoop(dp),
                            LatinSquareAdjacency, invt);
        end);
        
    # Latin square graphs from Cayley tables.
    InstallMethod(LatinSquareGraphCons,
        "for Cayley tables", true,
        [IsObject, IsList, IsBool], 1, function(filter, M, invt)
            return LatinSquareGraphCons(filter, QuasigroupByCayleyTable(M),
                                        invt);
        end);
fi;