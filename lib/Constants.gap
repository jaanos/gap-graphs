# The switching sets for Chang graphs
BindGlobal("ChangGraphSwitchingSet", [List([1..4], i -> [i, i+4]),
                            Set(List([1..8], i -> Set([i, (i mod 8)+1]))),
                            Union(List([1..3], i -> Set([i, (i mod 3)+1])),
                                List([1..5], i -> Set([i+3, (i mod 5)+4])))]);

# Involutions for the Chang graphs
BindGlobal("ChangGraphInvolution", [(), (1,8)(2,6)(3,7)(4,5), (1,3)(5,8)]);
