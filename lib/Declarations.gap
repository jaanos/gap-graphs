# Filters
DeclareFilter("NoVertexNames");
DeclareFilter("FullAutomorphismGroup");
DeclareFilter("IsSetGraph");
DeclareFilter("IsVectorGraph");
DeclareFilter("IsFormsGraph");
DeclareFilter("IsSpacesGraph");
DeclareFilter("IsConjugacyClassGraph");

# General constructions
DeclareConstructor("ProductGraphCons", [IsObject, IsList, IsFunction]);
DeclareConstructor("PowerGraphCons", [IsObject, IsRecord, IsInt, IsFunction]);
DeclareConstructor("GraphJoinCons", [IsObject, IsList]);
DeclareConstructor("ExtendedBipartiteDoubleGraphCons", [IsObject, IsRecord]);
DeclareConstructor("AntipodalQuotientGraphCons", [IsObject, IsRecord]);
DeclareConstructor("CliqueGraphCons", [IsObject, IsRecord, IsList]);
DeclareConstructor("IncidenceGraphCons", [IsObject, IsRecord, IsList]);

# Basic graphs
DeclareConstructor("CompleteMultipartiteGraphCons", [IsObject, IsInt, IsInt]);
DeclareConstructor("LatinSquareGraphCons", [IsObject, IsObject, IsBool]);
DeclareConstructor("HaarGraphCons", [IsObject, IsInt, IsList]);

# Set graphs
DeclareConstructor("KneserGraphCons", [IsObject, IsInt, IsInt, IsBool]);
DeclareConstructor("DoubledOddGraphCons", [IsObject, IsInt]);
DeclareConstructor("JohnsonGraphCons", [IsObject, IsInt, IsInt]);
DeclareConstructor("FoldedJohnsonGraphCons", [IsObject, IsInt]);
DeclareConstructor("ChangGraphCons", [IsObject, IsInt]);

# Vector graphs
DeclareConstructor("HammingGraphCons", [IsObject, IsInt, IsInt]);
DeclareConstructor("DoobGraphCons", [IsObject, IsInt, IsInt]);
DeclareConstructor("HalvedCubeGraphCons", [IsObject, IsInt]);
DeclareConstructor("FoldedCubeGraphCons", [IsObject, IsInt]);
DeclareConstructor("FoldedHalvedCubeGraphCons", [IsObject, IsInt]);
DeclareConstructor("BrouwerGraphCons", [IsObject, IsInt]);
DeclareConstructor("PasechnikGraphCons", [IsObject, IsInt]);
DeclareConstructor("AdditiveSymplecticCoverGraphCons",
                    [IsObject, IsInt, IsInt, IsInt]);
DeclareConstructor("MultiplicativeSymplecticCoverGraphCons",
                    [IsObject, IsInt, IsInt]);

# Forms graphs
DeclareConstructor("BilinearFormsGraphCons", [IsObject, IsInt, IsInt, IsInt]);
DeclareConstructor("HermiteanFormsGraphCons", [IsObject, IsInt, IsInt]);

# Spaces graphs
DeclareConstructor("GrassmannGraphCons", [IsObject, IsInt, IsInt, IsInt]);

# Named graphs
DeclareConstructor("ShrikhandeGraphCons", [IsObject]);
DeclareConstructor("PerkelGraphCons", [IsObject]);
DeclareConstructor("BiggsSmithGraphCons", [IsObject]);
