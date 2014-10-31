RequirePackage("grape");

CrossProduct := function(u, v)
    return [u[2]*v[3]-u[3]*v[2], u[3]*v[1]-u[1]*v[3], u[1]*v[2]-u[2]*v[1]];
end;

# The Brouwer graph Br(q) of pairs of 3-dimensional vectors over F_q,
# with two pairs being adjacent whenever the difference of the first vectors
# equals the cross product of the second vectors.
Brouwer := function(q)
    return Graph(Group(()), Elements(GF(q)^[2,3]), function(x,y) return x; end,
        function(x,y)
            return x <> y and x[1]-y[1] = CrossProduct(x[2], y[2]);
        end, true);
end;