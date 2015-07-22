# The tetrahedron with v=4, k=3, lm=2.
BindGlobal("TetrahedronGraph", CompleteGraph(SymmetricGroup(4)));

# The octahedron with v=6, k=4, lm=2, mu=4.
BindGlobal("OctahedronGraph", CocktailPartyGraph(3));

# The Petersen graph with v=10, k=3, lm=0, mu=1.
BindGlobal("PetersenGraph", OddGraph(2));

# The Clebsch graph with v=16, k=10, lm=6, mu=6.
BindGlobal("ClebschGraph", HalvedCubeGraph(5));

# The Schlaefli graph with v=27, k=16, lm=10, mu=8.
BindGlobal("SchlaefliGraph",
    Graph(DirectProduct(SymmetricGroup(6), SymmetricGroup(2)),
        [[-3,-3,1,1,1,1,1,1], [3,-1,-1,-1,-1,-1,-1,3]],
        OnRoots, RootAdjacency));

# The Hoffman-Singleton graph with v=50, k=7, lm=0, mu=1.
BindGlobal("HoffmanSingletonGraph", List([function()
        local G, dp, p1, p2, p3, p4;
        G := Group((1, 2, 3, 4, 5));
        dp := DirectProduct(G, G, G, Group((1, 2)));
        p1 := Projection(dp, 1);
        p2 := Projection(dp, 2);
        p3 := Projection(dp, 3);
        p4 := Projection(dp, 4);
        return Graph(dp, Cartesian(GF(5), GF(5), [1, 2]),
            function(x, g)
                local g4, j, s, u, v;
                u := 5^Image(p2, g) * Z(5)^0;
                v := 5^Image(p3, g) * Z(5)^0;
                g4 := Image(p4, g);
                s := x[3]^2;
                j := [1, 2*x[3]];
                return [1^g4*(x[1]+x[2]*(u-s*v)+s*(u^2-v^2)/2)+5^Image(p1, g),
                    j[1^g4]*(x[2]+s*u+v), x[3]^g4];
            end, function(x, y)
                return (x{[2,3]} = y{[2,3]} and y[1] in [x[1]-x[3], x[1]+x[3]])
                    or x[3] <> y[3] and y[1] = x[1] + x[2]*y[2]*x[3]^2;
            end, true);
    end])[1]());

# The Gewirtz graph with v=56, k=10, lm=0, mu=2.
BindGlobal("GewirtzGraph", Graph(MathieuGroup(21), [[1,2,3,7,10,20]],
                                    OnSets, DisjointSets));

# The strongly regular Witt graph with v=77, k=16, lm=0, mu=4.
BindGlobal("WittStronglyRegularGraph", Graph(MathieuGroup(22),
                                    [[1,2,3,7,10,20]], OnSets, DisjointSets));

# The graph with v=210, k=99, lm=48, mu=45 constructed by M. Klin
BindGlobal("KlinGraph", List([function()
        local G, H, N, V;
        H := Group([(1,2,3,4,5,6), (1,4)]);
        V := Set(List(SymmetricGroup(7), g -> H*g));
        N := Position(V, H*());
        G := EdgeOrbitsGraph(Action(SymmetricGroup(7), V, OnRight),
            List([(3,4,5,6,7), (2,4,6,3,5,7), (3,5,6,7), (2,4,5,6,7,3),
                        (2,3,5,6,7), (2,4,5,6), (2,3,5,6)],
                    g -> [N, Position(V, H*g)]));
        AssignVertexNames(G, V);
        return G;
    end])[1]());

# The Heawood graph with intersection array {3, 2, 2; 1, 1, 3}.
BindGlobal("HeawoodGraph", DesarguesianPlaneIncidenceGraph(2));

# The icosahedron with intersection array {5, 2, 1; 1, 2, 5}.
BindGlobal("IcosahedronGraph", MultiplicativeSymplecticCoverGraph(5, 2));

# The Sylvester graph with intersection array {5, 4, 2; 1, 1, 4}
BindGlobal("SylvesterGraph", Graph(SymmetricGroup(6),
                    [[1, [[[1,2], [3,4], [5,6]],
                          [[1,3], [2,5], [4,6]],
                          [[1,4], [2,6], [3,5]],
                          [[1,5], [2,4], [3,6]],
                          [[1,6], [2,3], [4,5]]]]],
                    function(x, g)
                        return [x[1]^g,
                            Set(List(x[2], l -> OnSetsSets(l, g)))];
                    end, function(x, y)
                        return x[1] <> y[1] and x[2] <> y[2] and
                            Set([x[1], y[1]]) in Intersection(x[2], y[2])[1];
                    end));

# The Perkel graph with intersection array {6, 5, 2; 1, 1, 3}.
BindGlobal("PerkelGraph", Graph(PSL(2, 19),
                    Elements(Filtered(ConjugacyClassesSubgroups(PSL(2, 19)),
                        x -> Size(x) = 57 and Order(x[1]) = 60)[1]),
                    ConjugateGroup, function(x, y)
                        return Order(Intersection(x, y)) = 10;
                    end, true));


# The Gosset graph with intersection array {27, 10, 1; 1, 10, 27}.
BindGlobal("GossetGraph",
    Graph(DirectProduct(SymmetricGroup(8), SymmetricGroup(2)),
        [[-3,-3,1,1,1,1,1,1]], OnRoots, RootAdjacency));

# The truncated Witt graph with intersection array {15, 14, 12; 1, 1, 9}.
BindGlobal("Witt23Graph", Graph(MathieuGroup(23), [[1,4,8,12,13,19,21,23]],
                                OnSets, DisjointSets));

# The large Witt graph with intersection array {30, 28, 24; 1, 3, 15}.
BindGlobal("Witt24Graph", Graph(MathieuGroup(24), [[1,4,8,12,13,19,21,23]],
                                OnSets, DisjointSets));

# The bipartite graph associated to Higman's design
# with intersection array {50, 49, 36; 1, 14, 50}.
BindGlobal("HigmanGraph", Graph(Stabilizer(MathieuGroup(24), [23, 24], OnSets),
                            [[1,4,8,12,13,19,21,23]], OnSets,
                            function (x, y)
                                return Length(Difference(Intersection(x, y),
                                                         [23, 24])) in [0, 4];
                            end));

# The Coxeter graph with intersection array {3,2,2,1; 1,1,1,2}.
BindGlobal("CoxeterGraph", PolarGraphNOorth(1, 3, 7));

# The doubly truncated Witt graph with intersection array {7,6,4,4; 1,1,1,6}.
BindGlobal("Witt22Graph", Graph(MathieuGroup(22), [[1,2,3,4,5,10,18,21]],
                                OnSets, DisjointSets));

# The dodecahedron with intersection array {3,2,1,1,1; 1,1,1,2,3}.
BindGlobal("DodecahedronGraph", List([function()
        local G, O, act, dp, p1, p2, u, v;
        dp := DirectProduct(AlternatingGroup(5), Group((1, 2)));
        p1 := Projection(dp, 1);
        p2 := Projection(dp, 2);
        act := function(t, g)
            return OnTuples(List([1,2], i -> t[i^Image(p2, g)]), Image(p1, g));
        end;
        u := [1, 2];
        v := [3, 4];
        O := Arrangements([1..5], 2);
        G := EdgeOrbitsGraph(Action(dp, O, act),
                             [Position(O, u), Position(O, v)]);
        AssignVertexNames(G, O);
        return G;
    end])[1]());

# The Desargues graph with intersection array {3,2,2,1,1; 1,1,2,2,3}.
BindGlobal("DesarguesGraph", DoubledOddGraph(2));

# The Biggs-Smith graph with intersection array {3,2,2,2,1,1,1; 1,1,1,1,1,1,3}.
BindGlobal("BiggsSmithGraph", Graph(PSL(2, 17),
                    Elements(Filtered(ConjugacyClassesSubgroups(PSL(2, 17)),
                        x -> Size(x) = 102 and Order(x[1]) = 24)[1]),
                    ConjugateGroup, function(x, y)
                        return Order(Intersection(x, y)) = 8;
                    end, true));
