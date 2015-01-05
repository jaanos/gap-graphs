# The Petersen graph.
BindGlobal("PetersenGraph", OddGraph(2));

# The Gewirtz graph.
BindGlobal("GewirtzGraph", Graph(MathieuGroup(21), [[1,2,3,7,10,20]],
                                    OnSets, DisjointSets));

# The Coxeter graph with intersection array {3,2,2,1; 1,1,1,2}
BindGlobal("CoxeterGraph", Graph(PGL(3, 2), [[1, 2, 4]],
                                    OnSets, DisjointSets));

# The Biggs-Smith graph with intersection array {3,2,2,2,1,1,1; 1,1,1,1,1,1,3}
BindGlobal("BiggsSmith", Graph(PSL(2, 17),
                    [Group([(1,16)(2,5)(4,13)(6,10)(7,12)(8,17)(9,14)(11,15),
                            (1,5,7,2)(3,10,4,15)(6,11,14,17)(8,18,9,13)])],
                    ConjugateGroup, function(x, y)
                        return Order(Intersection(x, y)) = 8;
                    end));
