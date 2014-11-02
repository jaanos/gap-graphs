# The Hermitean forms graph Her(d, r) of Hermitean matrices over GF(r^2).
BindGlobal("HermiteanFormsGraph", function(d, r)
    return Graph(Group(()), List(GF(r)^[d, d], x -> ToHermitean(x, r)),
        function(x, y) return x; end,
        function(x, y)
            return RankMat(x-y) = 1;
        end, true);
end);
