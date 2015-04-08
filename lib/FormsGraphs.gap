# The bilinear forms graph H_q(d, e) of matrices over GF(r^2).
BindGlobal("BilinearFormsGraph", function(q, d, e)
    local dp;
    dp := DirectProduct(Concatenation([GL(d, q), GL(e, q)],
        ListWithIdenticalEntries(d*e, FieldAdditionPermutationGroup(q))));
    return Graph(dp, Elements(GF(q)^[d,e]), OnMatrices(q, d, e, dp),
                RankAdjacency([1]), true);
end);

# The Hermitean forms graph Her(d, r) of Hermitean matrices over GF(r^2).
BindGlobal("HermiteanFormsGraph", function(d, r)
    local dp;
    dp := DirectProduct(Concatenation([GL(d, r^2)],
        ListWithIdenticalEntries(d*d, FieldAdditionPermutationGroup(r))));
    return Graph(dp, List(GF(r)^[d, d], x -> ToHermitean(x, r)),
        OnHermiteanMatrices(r, d, dp), RankAdjacency([1]), true);
end);
