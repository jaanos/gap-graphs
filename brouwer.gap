RequirePackage("grape");

Cross := function(u, v)
    return [u[2]*v[3]-u[3]*v[2], u[3]*v[1]-u[1]*v[3], u[1]*v[2]-u[2]*v[1]];
end;

Brouwer := function(q)
    return Graph(Group(()), Elements(GF(q)^[2,3]), function(x,y) return x; end,
        function(x,y)
            return x <> y and x[1]-y[1] = Cross(x[2], y[2]);
        end, true);
end;