# The bilinear forms graph H_q(d, e) of matrices over GF(r^2).
BindGlobal("BilinearFormsGraph", function(q, d, e)
    return Graph(WreathProduct(FieldAdditionPermutationGroup(q),
            MatrixColumnPermutationGroup(d, e)), Elements(GF(q)^[d,e]),
        OnMatrices(q, d, e), function(x, y)
            return RankMat(x-y) = 1;
        end, true);
end);

# The Hermitean forms graph Her(d, r) of Hermitean matrices over GF(r^2).
BindGlobal("HermiteanFormsGraph", function(d, r)
    return Graph(WreathProduct(FieldAdditionPermutationGroup(r),
            MatrixRowColumnPermutationGroup(d)),
        List(GF(r)^[d, d], x -> ToHermitean(x, r)),
        OnHermiteanMatrices(r, d), function(x, y)
            return RankMat(x-y) = 1;
        end, true);
end);
