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

# sparse6 string
BindGlobal("Sparse6String", function(G)
    local i, j, b, c, v, w;
    w := Log2Int(G.order-1)+1;
    b := [];
    v := 1;
    for i in [1..G.order] do
        for j in [1..i] do
            if IsVertexPairEdge(G, i, j) then
                if i = v+1 then
                    c := 1;
                else
                    if i > v then
                        Add(b, 1);
                        Append(b, IntToBits(i-1, w));
                    fi;
                    c := 0;
                fi;
                Add(b, c);
                Append(b, IntToBits(j-1, w));
                v := i;
            fi;
        od;
    od;
    if G.order in [2,4,8,16] and (-Length(b)) mod 6 > w
            and Length(Adjacency(G, G.order-1)) > 0
            and Length(Adjacency(G, G.order)) = 0 then
        Add(b, 0);
    fi;
    return List(Concatenation([58], Graph6EncodeInteger(G.order),
                                BitsToString(b, 1)), CharInt);
end);

# Read graph from sparse6 string
BindGlobal("GraphFromSparse6String", function(s)
    local A, b, t, i, m, n, v, w, c, x, y, z;
    if s[1] = ':' then
        s := s{[2..Length(s)]};
    fi;
    s := List(s, IntChar);
    t := Graph6DecodeInteger(1, s);
    i := t[1];
    n := t[2];
    w := Log2Int(n-1)+1;
    b := StringToBits(s{[i..Length(s)]});
    A := List([1..n], j -> []);
    i := 1;
    v := 1;
    z := BitsToInt(b, w+1);
    m := 2^w;
    for y in z do
        c := Int(y/m);
        x := y mod m + 1;
        if c = 1 then
            v := v+1;
        fi;
        if v > n or x > n then
            break;
        fi;
        if x > v then
            v := x;
        else
            Add(A[x], v);
            Add(A[v], x);
        fi;
    od;
    return AdjFunGraph([1..n], ListAdjacency(A));
end);

# sparse6 string
BindGlobal("IncrementalSparse6String", function(G, H)
    local i, j, b, c, v, w;
    if G.order <> H.order then
        Error("graphs have different orders");
        return fail;
    fi;
    w := Log2Int(G.order-1)+1;
    b := [];
    v := 1;
    for i in [1..G.order] do
        for j in [1..i] do
            if IsVertexPairEdge(G, i, j) <> IsVertexPairEdge(H, i, j) then
                if i = v+1 then
                    c := 1;
                else
                    if i > v then
                        Add(b, 1);
                        Append(b, IntToBits(i-1, w));
                    fi;
                    c := 0;
                fi;
                Add(b, c);
                Append(b, IntToBits(j-1, w));
                v := i;
            fi;
        od;
    od;
    if G.order in [2,4,8,16] and (-Length(b)) mod 6 > w
            and Adjacency(G, G.order-1) = Adjacency(H, G.order-1)
            and Adjacency(G, G.order) = Adjacency(H, G.order) then
        Add(b, 0);
    fi;
    return List(Concatenation([59], BitsToString(b, 1)), CharInt);
end);

# Read graph from incremental sparse6 string
BindGlobal("GraphFromIncrementalSparse6String", function(s, H)
    local A, B, b, i, m, v, w, c, x, y, z;
    if s[1] = ';' then
        s := s{[2..Length(s)]};
    fi;
    B := List([1..H.order], j -> Adjacency(H, j));
    A := List(B, ShallowCopy);
    s := List(s, IntChar);
    w := Log2Int(H.order-1)+1;
    b := StringToBits(s);
    i := 1;
    v := 1;
    z := BitsToInt(b, w+1);
    m := 2^w;
    for y in z do
        c := Int(y/m);
        x := y mod m + 1;
        if c = 1 then
            v := v+1;
        fi;
        if v > H.order or x > H.order then
            break;
        fi;
        if x > v then
            v := x;
        elif x in A[v] then
            Remove(A[x], Position(A[x], v));
            if x <> v then
                Remove(A[v], Position(A[v], x));
            fi;
        else
            Add(A[x], v);
            Add(A[v], x);
        fi;
    od;
    return AdjFunGraph([1..H.order], ListAdjacency(A));
end);

# auto6 string
BindGlobal("Auto6String", function(G)
    local b, g, h, i, l, p, r, s, w, x, sch;
    g := GeneratorsOfGroup(G.group);
    l := [1..Length(g)];
    i := 0;
    for p in [1..Length(g)] do
        if g[p] = () then
            i := i+1;
        else
            l[p] := l[p] - i;
        fi;
    od;
    sch := List(G.schreierVector, function(y)
                                    if y < 0 then
                                        return y;
                                    else
                                        return l[y];
                                    fi;
                                  end);
    g := Filtered(g, y -> y <> ());
    w := Log2Int(G.order-1)+1;
    r := Length(G.representatives);
    b := IntToBits(r, w);
    for i in [1..r] do
        Append(b, IntToBits(G.representatives[i]-1, w));
        Append(b, IntToBits(Length(G.adjacencies[i]), w));
        for x in G.adjacencies[i] do
            Append(b, IntToBits(x-1, w));
        od;
    od;
    for h in g do
        for i in [1..G.order] do
            Append(b, IntToBits(i^h - 1, w));
        od;
    od;
    if Length(g) > 1 then
        s := Log2Int(Length(g)) + 1;
        for i in [1..G.order] do
            if sch[i] < 0 then
                Append(b, IntToBits(0, s));
            else
                Append(b, IntToBits(sch[i], s));
            fi;
        od;
    fi;
    return List(Concatenation([33], Graph6EncodeInteger(G.order),
                Graph6EncodeInteger(Length(g)), BitsToString(b, 0)), CharInt);
end);

# Read graph from auto6 string
BindGlobal("GraphFromAuto6String", function(s)
    local G, A, R, S, a, b, d, r, t, g, i, j, k, n, v, w;
    if s[1] = '!' then
        s := s{[2..Length(s)]};
    fi;
    s := List(s, IntChar);
    t := Graph6DecodeInteger(1, s);
    i := t[1];
    n := t[2];
    w := Log2Int(n-1)+1;
    t := Graph6DecodeInteger(i, s);
    i := t[1];
    g := t[2];
    b := StringToBits(s{[i..Length(s)]});
    d := BitsToInt(b, w);
    r := d[1];
    if r = 0 then
        r := n;
    fi;
    i := 2;
    R := [];
    A := [];
    for j in [1..r] do
        Add(R, d[i]+1);
        Add(A, d{[i+2..i+d[i+1]+1]}+1);
        i := i + d[i+1] + 2;
    od;
    if g = 0 then
        G := ();
        v := List([1..n], x -> -Position(R, x));
    else
        G := [];
        for j in [1..g] do
            Add(G, PermList(d{[i..i+n-1]}+1));
            i := i+n;
        od;
        if g = 1 then
            v := List([1..n], function(x)
                local y;
                y := Position(R, x);
                if y = fail then
                    return 1;
                else
                    return -y;
                fi;
            end);
        else
            k := w*(i-1);
            t := Log2Int(g)+1;
            v := BitsToInt(b{[k+1..k+n*t]}, t);
            for j in [1..r] do
                v[R[j]] := -j;
            od;
        fi;
    fi;
    return rec(
        isGraph := true,
        order := n,
        names := [1..n],
        representatives := Immutable(R),
        adjacencies := A,
        group := Group(G),
        schreierVector := Immutable(v)
    );
end);
