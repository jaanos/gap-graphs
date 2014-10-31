RequirePackage("grape");

V := GF(4)^4;
P := Elements(Subspaces(V, 1));
LP := Elements(Subspaces(V, 2));
LLP := List(LP, l -> Filtered([1..Size(P)], i -> IsSubset(l, P[i])));
O := [1, 2, 6, 11, 17, 20];
GQ35 := Graph(Group(()), [22..85], function(x,y) return x; end,
    function(x, y)
        local l, ll;
        if x = y then
            return false;
        fi;
        ll := Filtered(LLP, l -> IsSubset(l, [x, y]));
        return Size(ll) = 1 and Size(Intersection(ll[1], O)) = 1;
    end, true);;
    
spr35 := [ [ 1, 2, 3, 4 ], [ 5, 6, 7, 8 ], [ 9, 10, 11, 12 ], [ 13, 14, 15, 16 ], [ 17, 18, 19, 20 ], [ 21, 22, 23, 24 ], [ 25, 26, 27, 28 ], [ 29, 30, 31, 32 ], 
  [ 33, 34, 35, 36 ], [ 37, 38, 39, 40 ], [ 41, 42, 43, 44 ], [ 45, 46, 47, 48 ], [ 49, 50, 51, 52 ], [ 53, 54, 55, 56 ], [ 57, 58, 59, 60 ], 
  [ 61, 62, 63, 64 ] ];;
  
GQ35spr := Graph(Group(()), [1..64], function(x,y) return x; end,
    function(x,y)
        return Distance(GQ35, x, y) = 1 and Size(Filtered(spr35, s -> x in s and y in s)) = 0;
    end, true);;
    
GQ53 := Graph(Group(()), Cliques(GQ35, 4), function(x,y) return x; end,
    function(x,y)
        return Size(Intersection(x, y)) = 1;
    end, true);;
    
spr53 := [ [ 1, 2, 3, 4, 5, 6 ], [ 7, 22, 27, 28, 29, 30 ], [ 8, 33, 55, 58, 81, 83 ], [ 9, 34, 54, 57, 85, 87 ], [ 10, 36, 43, 65, 93, 94 ], 
  [ 11, 35, 42, 66, 90, 95 ], [ 12, 39, 48, 49, 50, 51 ], [ 13, 26, 53, 63, 89, 92 ], [ 14, 25, 52, 64, 91, 96 ], [ 15, 37, 46, 59, 75, 79 ], 
  [ 16, 38, 47, 60, 74, 80 ], [ 17, 56, 69, 70, 71, 72 ], [ 18, 31, 40, 61, 73, 77 ], [ 19, 32, 41, 62, 76, 78 ], [ 20, 23, 45, 68, 84, 88 ], 
  [ 21, 24, 44, 67, 82, 86 ] ];;
  
GQ53spr := Graph(Group(()), [1..96], function(x,y) return x; end,
    function(x,y)
        return Distance(GQ53, x, y) = 1 and Size(Filtered(spr53, s -> x in s and y in s)) = 0;
    end, true);;