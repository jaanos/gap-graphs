# Filters
DeclareFilter("NoVertexNames");
DeclareFilter("FullAutomorphismGroup");
DeclareFilter("IsSetGraph");
DeclareFilter("IsVectorGraph");
DeclareFilter("IsConjugacyClassGraph");

# General constructions
DeclareConstructor("ProductGraphCons", [IsObject, IsList, IsFunction]);
DeclareConstructor("GraphJoinCons", [IsObject, IsList]);
DeclareConstructor("ExtendedBipartiteDoubleGraphCons", [IsObject, IsRecord]);
DeclareConstructor("AntipodalQuotientGraphCons", [IsObject, IsRecord]);
DeclareConstructor("CliqueGraphCons", [IsObject, IsRecord, IsList]);
DeclareConstructor("IncidenceGraphCons", [IsObject, IsRecord, IsList]);

# Basic graphs
DeclareConstructor("CompleteMultipartiteGraphCons", [IsObject, IsInt, IsInt]);
DeclareConstructor("LatinSquareGraphCons", [IsObject, IsObject, IsBool]);

# Set graphs
DeclareConstructor("JohnsonGraphCons", [IsObject, IsInt, IsInt]);
DeclareConstructor("FoldedJohnsonGraphCons", [IsObject, IsInt]);

# Vector graphs
DeclareConstructor("HammingGraphCons", [IsObject, IsInt, IsInt]);

# Named graphs
DeclareConstructor("ShrikhandeGraphCons", [IsObject]);
DeclareConstructor("PerkelGraphCons", [IsObject]);
DeclareConstructor("BiggsSmithGraphCons", [IsObject]);
