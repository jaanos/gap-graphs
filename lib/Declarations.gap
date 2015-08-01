# Filters
DeclareFilter("NoVertexNames");
DeclareFilter("FullAutomorphismGroup");
DeclareFilter("IsVectorGraph");

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

# Vector graphs
DeclareConstructor("HammingGraphCons", [IsObject, IsInt, IsInt]);
DeclareConstructor("ShrikhandeGraphCons", [IsObject]);