# The auto6 graph format

**auto6** is a format for storing graphs without multiple edges in a compact manner, using only printable ASCII characters. Files in this format have text type and contain one line per graph.

**auto6** has been inspired by Brendan McKay's [**graph6** and **sparse6**](https://cs.anu.edu.au/~bdm/data/formats.html) formats, and also by GRAPE's representation of graphs using permutation groups. As such, **auto6** is suitable for graphs with many symmetries, i.e. a small number of vertex orbits.

## General principles

As in the **graph6** and **sparse6** formats, there is, apart from the header, one object per line. Apart from  the header, end-of-line characters, and the character `!` which starts a line, all bytes have a value in the range 63-126 (which are all printable ASCII characters). A file of objects is a text file, so whatever end-of-line convention is locally used is fine.

In the following definitions the functions `R(x)` (byte representation of a bitstring `x`) and `N(n)` (byte representation of an integer n in the range 0-68719476735 (2<sup>36</sup>-1)) are the same as those used in [**graph6** and **sparse6**](https://cs.anu.edu.au/~bdm/data/formats.txt).

Besides the vertices and the edges of the graph, **auto6** stores two additional pieces of information: a group of automorphisms of the graph (which need not be the full group of automorphisms), and the [Schreier vector](https://en.wikipedia.org/wiki/Schreier_vector), which tells which of the generators has been used to reach a given vertex from the orbit representatives.

## Description of the auto6 format

* Data type: graphs of order 0 to 68719476735. Directed edges and loops are permitted, while multiple edges are not.

* Optional header:
```
>>auto6<<
```
(without end of line!)
    
* File name extension: `.a6`

* General structure:

Each graph occupies one text line. Except for the first character and end-of-line characters, each byte has the form 63+x, where 0 <= x <= 63. The byte encodes the six bits of x.

The encoded graph consists of:

1. The character `!`. (This is present to distinguish the code from **graph6** and **sparse6** formats.)

2. The number of vertices.

3. The number of generators of a group.

4. The number of vertex orbit representatives.

5. For each vertex orbit representative:
  - its number,
  - the number of (out-)neighbours,
  - a list of neighbours.
  
6. A list of group generators.

7. The Schreier vector.

8. End-of-line.

### Number of vertices `n`

1, 4, or 8 bytes `N(n)`.
This is the same as **graph6** and **sparse6** format.

### Number of generators `g`

1, 4, or 8 bytes `N(g)`.
If the group is trivial, then `g` is zero (note that in this case the format is a less efficient version of **sparse6**).


The remaining bytes are `R(x)`, where `x` is a bitstring obtained by concatenating the bitstrings described below. `k` will denote the number of bits needed to represent `n-1` in binary.

### Number of representatives `r`

`k` bits. When `n` is a power of 2 and the group is trivial (i.e., `g` is zero), then we have `r = n`, which is encoded as `k` zero bits. Otherwise, encode `r` in bigendian order.

### Vertex orbit representative `i`

`k` bits. Encode the index `i` of the representative (an integer between 0 and `n-1`) in bigendian order.

#### Number of neighbours `a[i]`

`k` bits. Encode the number `a[i]` in bigendian order.

#### List of neighbours

`a[i]*k` bits. Encode the index of each neighbour in bigendian order.

### List of generators

`g*n*k` bits. Write each generator as a permutation of vertex indices (integers between 0 and `n-1`); encode each index in bigendian order.

### Schreier vector

Let `t` be the number of bits needed to represent `g` in binary.
`n*t` bits. Encode each entry in the Schreier vector in bigendian order, substituting negative entries (at representatives) by zeros.

If the group is trivial (i.e., `g` is zero), the Schreier vector is omitted.

## Example

    !IACBHFCcTfAHKBGSaVT`XTg
    
`!` indicates **auto6** format.
Subtract 63 from the other bytes and write them in binary,  six bits each.

    001010 000010 000100 000011 001001 000111 000100 100100
    010101 100111 000010 001001 001100 000011 001000 010100
    100010 010111 010101 100001 011001 010101 101000
    
The first byte is not 63, so it is `n = 10`.
The second byte is not 63, so it is `g = 2`.
`n-1` needs `k = 4` bits. Write the following bits in groups of `k`:

    0001 0000 0011 0010 0100 0111
    0001 0010 0100 0101 0110 0111 0000 1000 1001 0011
    0000 0011 0010 0001 0100 1000 1001 0111 0101 0110
    
The first nibble is the number of representatives `r = 1`.
The only representative has index 0 and 3 neighbours, namely 2, 4, and 7.
The next two lines are the two generators, namely the permutations

    (0 1 2 4 6)(3 5 7 8 9)
    (1 3)(5 8)(6 9)
    
The two permutations generate a group isomorphic to the symmetric group on 5 elements. The graph being defined here is the Petersen graph.

Finally, `g` needs `t = 2` bits, so write the remaining bits in groups of `t`:

    00 01 01 10 01 01 01 01 10 10 00

The Schreier vector with `n` elements is

    [-1, 1, 1, 2, 1, 1, 1, 1, 2, 2]
    
The last two bits are just padding.
