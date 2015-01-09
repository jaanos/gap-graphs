# The Petersen graph with v=10, k=3, lm=0, mu=1.
BindGlobal("PetersenGraph", OddGraph(2));

# The Clebsch graph with v=16, k=10, lm=6, mu=6.
BindGlobal("ClebschGraph", HalvedCubeGraph(5));

# The Schlaefli graph with v=27, k=16, lm=10, mu=8.
BindGlobal("SchlaefliGraph",
    Graph(DirectProduct(SymmetricGroup(6), SymmetricGroup(2)),
        [[-3,-3,1,1,1,1,1,1], [3,-1,-1,-1,-1,-1,-1,3]],
        OnRoots, RootAdjacency));

# The Gewirtz graph with v=56, k=10, lm=0, mu=2.
BindGlobal("GewirtzGraph", Graph(MathieuGroup(21), [[1,2,3,7,10,20]],
                                    OnSets, DisjointSets));

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

# The icosahedron.
BindGlobal("IcosahedronGraph", MultiplicativeSymplecticCoverGraph(5, 2));

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

# The Coxeter graph with intersection array {3,2,2,1; 1,1,1,2}.
BindGlobal("CoxeterGraph", Graph(PGL(3, 2), [[1, 2, 4]],
                                    OnSets, DisjointSets));

# The Biggs-Smith graph with intersection array {3,2,2,2,1,1,1; 1,1,1,1,1,1,3}.
BindGlobal("BiggsSmithGraph", Graph(PSL(2, 17),
                    Elements(Filtered(ConjugacyClassesSubgroups(PSL(2, 17)),
                        x -> Size(x) = 102 and Order(x[1]) = 24)[1]),
                    ConjugateGroup, function(x, y)
                        return Order(Intersection(x, y)) = 8;
                    end, true));
