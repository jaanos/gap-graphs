gap-graphs
==========

Various constructions of graphs that I have encountered in my research.

For now, these are standalone GAP scripts requiring the GRAPE package. Eventually, I intend to reorganize them into a proper GAP package.

The following graphs are currently supported:

* Classical distance-regular graphs:
  * Dual polar graphs B<sub>d</sub>(q)
  * Dual polar graphs C<sub>d</sub>(q)
  * Dual polar graphs D<sub>d</sub>(q)
* Other infinite families of distance-regular graphs:
  * the Brouwer "cross product" graphs Br(q) with intersection arrays {q<sup>3</sup>-1, q<sup>3</sup>-q, q<sup>3</sup>-q<sup>2</sup>+1; 1, q, q<sup>2</sup>-1}
  * de Caen, Mathon and Moorhouse's Preparata graphs Pr(t, h, e) with intersection arrays {2<sup>2t</sup>-1, 2<sup>2t</sup>-2<sup>h</sup>, 1; 1, 2<sup>h</sup>, 2<sup>2t</sup>-1}
* Sporadic distance-regular graphs:
  * the Doro graph on 65 vertices with intersection array {10, 6, 4; 1, 2, 5}
  * the Doro graph on 68 vertices with intersection array {12, 10, 3; 1, 3, 8}
  * the Coxeter graph with intersection array {3, 2, 2, 1; 1, 1, 1, 2}
  * the Biggs-Smith graph with intersection array {3, 2, 2, 2, 1, 1, 1; 1, 1, 1, 1, 1, 1, 3}
* Strongly regular graphs:
  * the Gewirtz graph with v = 56, k = 10, λ = 0, μ = 2
