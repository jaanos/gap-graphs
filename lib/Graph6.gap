# Constants
BindGlobal("Graph6ByteInteger", 62);
BindGlobal("Graph6WordInteger", 258047);
BindGlobal("Graph6DwordInteger", 68719476735);
BindGlobal("Graph6MaxInteger", Graph6DwordInteger);
BindGlobal("Graph6BitVector", List([0..5], i -> 2^(5-i)));
BindGlobal("Graph6IntVector", List([0..5], i -> 64^(5-i)));

# graph6 integer encoding
BindGlobal("Graph6EncodeInteger", function(x)
    if not IsInt(x) or x < 0 or x >= Graph6MaxInteger then
        Error("not an integer in allowed range");
        return fail;
    fi;
    if x <= Graph6ByteInteger then
        return [x+63];
    elif x <= Graph6WordInteger then
        return [63, Int(x/4096), Int(x/64), x] mod 64 + 63;
    else
        return Concatenation([126, 126],
                List(Graph6IntVector, y -> Int(x/y) mod 64 + 63));
    fi;
end);

# graph6 integer decoding
BindGlobal("Graph6DecodeInteger", function(i, s)
    if s[i] <> 126 then
        return [i+1, s[i] - 63];
    elif s[i+1] <> 126 then
        return [i+4, (s{[i+1..i+3]} - 63)*Graph6IntVector{[4..6]}];
    else
        return [i+8, (s{[i+2..i+7]} - 63)*Graph6IntVector];
    fi;
end);

# Convert an integer to a bitstring of length w
BindGlobal("IntToBits", function(x, w)
    return List([1..w], i -> Int(x / 2^(w-i)) mod 2);
end);

# Convert a bitstring to a list of integers given the width w
BindGlobal("BitsToInt", function(b, w)
    local out, i;
    out := [];
    i := 0;
    while i+w <= Length(b) do
        Add(out, b{[i+1..i+w]} * List([1..w], j -> 2^(w-j)));
        i := i+w;
    od;
    return out;
end);

# Convert a bitstring to a graph6 string with padding p
BindGlobal("BitsToString", function(b, p)
    local c;
    c := Concatenation(b, ListWithIdenticalEntries(-Length(b) mod 6, p));
    return List([0,6..Length(c)-6], i -> c{[i+1..i+6]}*Graph6BitVector + 63);
end);

# Convert a graph6 string to a bitstring
BindGlobal("StringToBits", function(s)
    return Concatenation(List(s, x -> List(Graph6BitVector, y -> Int((x-63)/y) mod 2)));
end);

# graph6 string
BindGlobal("Graph6String", function(G)
    local A, s;
    A := CollapsedAdjacencyMat(Group(()), G);
    s := Concatenation(Graph6EncodeInteger(G.order),
        BitsToString(Concatenation(List([1..G.order], i -> A[i]{[1..i-1]})), 0));
    return List(s, CharInt);
end);

# Read graph from graph6 string
BindGlobal("GraphFromGraph6String", function(s)
    local A, b, n, i, j, t;
    s := List(s, IntChar);
    t := Graph6DecodeInteger(1, s);
    i := t[1];
    n := t[2];
    b := StringToBits(s{[i..Length(s)]});
    A := NullMat(n, n);
    i := 1;
    for j in [2..n] do
        A[j]{[1..j-1]} := b{[i..i+j-2]};
        A{[1..j-1]}[j] := b{[i..i+j-2]};
        i := i+j-1;
    od;
    return AdjFunGraph([1..n], MatrixAdjacency(A));
end);
