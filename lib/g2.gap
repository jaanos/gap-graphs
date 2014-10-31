RequirePackage("grape");

GC := function(x,y,i,j) return x[i]*y[j] - x[j]*y[i]; end;

G2 := function(q)
    local P, Q;
    if q mod 2 = 0 then
        Q := [[0,0,0,0,1,0,0],[0,0,0,0,0,1,0],[0,0,0,0,0,0,1],[0,0,0,1,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]]*Z(q)^0;;
    else
        Q := [[0,0,0,0,1,0,0],[0,0,0,0,0,1,0],[0,0,0,0,0,0,1],[0,0,0,1,0,0,0],[1,0,0,0,0,0,0],[0,1,0,0,0,0,0],[0,0,1,0,0,0,0]]*Z(q)^0;;
    fi;
    P := Filtered(List(Subspaces(GF(q)^7, 1), y -> Elements(y)[2]), x -> x*Q*x = 0*Z(q));;

    return Graph(Group(()), P, function(x,y) return x; end,
        function(x,y)
            return x <> y and x*Q*y = 0*Z(q) and y*Q*x = 0*Z(q) and 
            [ GC(x,y,2,3), GC(x,y,6,5), GC(x,y,3,1), GC(x,y,7,6), GC(x,y,1,2), GC(x,y,5,7) ] =
            [ GC(x,y,4,5), GC(x,y,4,3), GC(x,y,4,6), GC(x,y,4,1), GC(x,y,4,7), GC(x,y,4,2) ];
        end, true);;
end;

