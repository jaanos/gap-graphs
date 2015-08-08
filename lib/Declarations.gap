# Filters
DeclareFilter("NoVertexNames");
DeclareFilter("FullAutomorphismGroup");
DeclareFilter("IsSetGraph");
DeclareFilter("IsVectorGraph");
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

# Set graphs
DeclareConstructor("KneserGraphCons", [IsObject, IsInt, IsInt, IsBool]);
DeclareConstructor("DoubledOddGraphCons", [IsObject, IsInt]);
DeclareConstructor("JohnsonGraphCons", [IsObject, IsInt, IsInt]);
DeclareConstructor("FoldedJohnsonGraphCons", [IsObject, IsInt]);
DeclareConstructor("ChangGraphCons", [IsObject, IsInt]);

# Vector graphs
DeclareConstructor("HammingGraphCons", [IsObject, IsInt, IsInt]);
DeclareConstructor("DoobGraphCons", [IsObject, IsInt, IsInt]);
DeclareConstructor("BrouwerGraphCons", [IsObject, IsInt]);
DeclareConstructor("PasechnikGraphCons", [IsObject, IsInt]);

# Named graphs
DeclareConstructor("ShrikhandeGraphCons", [IsObject]);
DeclareConstructor("PerkelGraphCons", [IsObject]);
DeclareConstructor("BiggsSmithGraphCons", [IsObject]);
