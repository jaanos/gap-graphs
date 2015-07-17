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
        return [63, 63, Int(x/1073741824), Int(x/16777216), Int(x/262144),
                Int(x/4096), Int(x/64), x] mod 64 + 63;
    fi;
end);

# graph6 integer decoding
BindGlobal("Graph6DecodeInteger", function(b, i)
    if b[i] <> 126 then
        return [b[i] - 63, i+1];
    elif b[i+1] <> 126 then
        return [(b{[i+1..i+3]} - 63)*Graph6IntVector{[4..6]}, i+4];
    else
        return [(b{[i+2..i+7]} - 63)*Graph6IntVector, i+8];
    fi;
end);
