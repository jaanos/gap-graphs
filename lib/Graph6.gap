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
BindGlobal("Graph6DecodeInteger", function(i, b)
    if b[i] <> 126 then
        return [i+1, b[i] - 63];
    elif b[i+1] <> 126 then
        return [i+4, (b{[i+1..i+3]} - 63)*Graph6IntVector{[4..6]}];
    else
        return [i+8, (b{[i+2..i+7]} - 63)*Graph6IntVector];
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
