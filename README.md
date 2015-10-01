gap-graphs
==========

Various constructions of graphs that I have encountered in my research.

Currently, the scripts are undergoing reorganization into a proper GAP package.
The GRAPE package for GAP is required.

The following graphs are currently supported:

* Classical distance-regular graphs:
  * Hamming graphs H(d, e) with hypercubes H(d, 2) as special cases
  * Halved cubes 1/2 H(d, 2)
  * Doob graphs Doob(n, d) of diameter 2n+d
  * Bilinear forms graphs H<sub>q</sub>(d, e)
  * Hermitean forms graphs Her(d, r)
  * Grassmann graphs J<sub>q</sub>(n, d)
  * Twisted Grassmann graphs TG<sub>d</sub>(q)
  * Dual polar graphs B<sub>d</sub>(q)
  * Dual polar graphs C<sub>d</sub>(q)
  * Dual polar graphs D<sub>d</sub>(q)
  * Dual polar graphs <sup>2</sup>D<sub>d+1</sub>(q)
  * Dual polar graphs <sup>2</sup>A<sub>e</sub>(q) of diameter [e/2]
  * the Gosset graph with intersection array {27, 10, 1; 1, 10, 27}
  * the large Witt graph with intersection array {30, 28, 24; 1, 3, 15}
  * the truncated Witt graph with intersection array {15, 14, 12; 1, 1, 9}
* Distance-regular graphs associated to finite geometries:
  * incidence graphs of Desarguesian projective planes
  * incidence graphs of Hall projective planes
  * incidence graphs of Hughes projective planes (both ordinary and exceptional)
  * collinearity graphs of generalized quadrangles Q(d, q) of order (q, q<sup>d-3</sup>) for d = 3, 4, 5
  * collinearity graphs of generalized quadrangles H(d, r) of order (r<sup>2</sup>, r<sup>d-5/2</sup>) for d = 3, 4
  * collinearity graphs of generalized quadrangles W(q) of order (q, q)
  * collinearity graphs of generalized quadrangles AS(q) of order (q-1, q+1)
  * collinearity graphs of derived generalized quadrangles T<sub>d</sub>(O) of order (q, q<sup>d-1</sup>)
  * collinearity graphs of derived generalized quadrangles T<sup>*</sup>(O) of order (2<sup>h</sup>-1, 2<sup>h</sup>+1)
  * collinearity graphs of derived generalized quadrangles P(S, x) of order (s-1, s+1)
* Other infinite families of distance-regular graphs:
  * Cycles C(n)
  * Complete Taylor graphs, i.e. K<sub>n,n</sub> minus a matching
  * Odd graphs O(d)
  * Folded Johnson graphs J̃(2d, d)
  * Folded cubes H̃(d, 2)
  * Folded halved cubes 1/2 H̃(2d, 2)
  * Doubled Odd graphs 2J(2d+1, d)
  * Doubled Grassmann graphs 2J<sub>q</sub>(2d+1, d)
  * unitary nonisotropics graphs with intersection arrays {q<sup>2</sup>-q, q<sup>2</sup>-q-2, q+1; 1, 1, q<sup>2</sup>-2q}
  * the Brouwer "vector product" graphs Br(q) with intersection arrays {q<sup>3</sup>-1, q<sup>3</sup>-q, q<sup>3</sup>-q<sup>2</sup>+1; 1, q, q<sup>2</sup>-1}
  * the Pasechnik graphs Pa(q) with intersection arrays {q<sup>3</sup>, q<sup>3</sup>-1, q<sup>3</sup>-q, q<sup>3</sup>-q<sup>2</sup>+1; 1, q, q<sup>2</sup>-1, q<sup>3</sup>}
  * coset graphs of Kasami codes
  * de Caen, Mathon and Moorhouse's Preparata graphs Pr(t, e) (h = 1) and their quotients (h > 1) with intersection arrays {2<sup>2t</sup>-1, 2<sup>2t</sup>-2<sup>h</sup>, 1; 1, 2<sup>h</sup>, 2<sup>2t</sup>-1}
  * additive symplectic covers of complete graphs (j = i) and their quotients (j < i) with intersection arrays {p<sup>ni</sup>-1, p<sup>ni</sup>-p<sup>ni-j</sup>, 1; 1, p<sup>ni-j</sup>, p<sup>ni</sup>-1}
  * multiplicative symplectic covers of complete graphs with intersection arrays {q, q-m-1, 1; 1, m, q} (q or m even, m divides q-1)
* Families of strongly regular graphs:
  * Cocktail party graphs
  * Paley graphs
  * Latin square graphs
  * polar graphs O<sup>(±)</sup>(d, q)
  * polar graphs NO<sup>±</sup><sub>⊥</sub>(d, q)
  * polar graphs Sp(d, q)
  * polar graphs U(d, r)
  * affine polar graphs VO<sup>±</sup>(2e, q)
  * affine polar graphs VNO<sup>±</sup>(2e, q)
* Sporadic and other named distance-regular graphs:
  * the Heawood graph with intersection array {3, 2, 2; 1, 1, 3}
  * the skeleton of the cube with intersection array {3, 2, 1; 1, 2, 3}
  * the skeleton of the icosahedron with intersection array {5, 2, 1; 1, 2, 5}
  * the Sylvester graph with intersection array {5, 4, 2; 1, 1, 4}
  * the Perkel graph with intersection array {6, 5, 2; 1, 1, 3}
  * the Doro graph on 65 vertices with intersection array {10, 6, 4; 1, 2, 5}
  * the Doro graph on 68 vertices with intersection array {12, 10, 3; 1, 3, 8}
  * the bipartite graph associated to Higman's design with intersection array {50, 49, 36; 1, 14, 50}
  * the Coxeter graph with intersection array {3, 2, 2, 1; 1, 1, 1, 2}
  * the doubly truncated Witt graph with intersection array {7, 6, 4, 4; 1, 1, 1, 6}
  * the skeleton of the dodecahedron with intersection array {3, 2, 1, 1, 1; 1, 1, 1, 2, 3}
  * the Desargues graph with intersection array {3, 2, 2, 1, 1; 1, 1, 2, 2, 3}
  * the Biggs-Smith graph with intersection array {3, 2, 2, 2, 1, 1, 1; 1, 1, 1, 1, 1, 1, 3}
* Named strongly regular graphs:
  * the skeleton of the tetrahedron with v = 4, k = 3, λ = 2
  * the skeleton of the octahedron with v = 6, k = 4, λ = 2, μ = 4
  * the Petersen graph with v = 10, k = 3, λ = 0, μ = 1
  * the Shrikhande graph with v = 16, k = 6, λ = 2, μ = 2
  * the Clebsch graph with v = 16, k = 10, λ = 6, μ = 6
  * the Schläfli graph with v = 27, k = 16, λ = 10, μ = 8
  * the three Chang graphs with v = 28, k = 12, λ = 6, μ = 4
  * the Hoffman-Singleton graph with v = 50, k = 7, λ = 0, μ = 1
  * the Gewirtz graph with v = 56, k = 10, λ = 0, μ = 2
  * the graph with v = 77, k = 16, λ = 0, μ = 4 associated to S(3, 6, 22)
  * the graph with v = 210, k = 99, λ = 48, μ = 45 constructed by M. Klin
* Other graphs:
  * Complete multipartite graphs
  * Haar graphs
  * Kneser graphs Kn(n, k)
* Graph derivation:
  * Box (Cartesian) product and power
  * Cross (tensor) product and power
  * Strong product and power
  * Graph join
  * Bipartite double
  * Extended bipartite double
  * Halved graphs
  * Antipodal quotients
  * Clique graphs
  * Vertex-clique incidence graphs
* Graphs from files:
  * incidence graphs from files as listed on E. Moorhouse's webpage on [projective planes](http://www.uwyo.edu/moorhouse/pub/planes/) or [generalized polygons](http://www.uwyo.edu/moorhouse/pub/genpoly/)
  * collinearity graphs from files as listed on E. Moorhouse's webpage on [generalized polygons](http://www.uwyo.edu/moorhouse/pub/genpoly/)
  * graphs in [graph6, sparse6](https://cs.anu.edu.au/~bdm/data/formats.html) and [auto6](auto6.md) formats

Some constructions, such as complete graphs, Johnson graphs and Cayley graphs,
are already available in GRAPE.
