////////////////
//  X_1(2,18) //
////////////////

load "coleman.m";
Q:=9*x^6*y^3 + 12*x^4*y^4 - 9*x^4*y^3 - 81*x^4*y^2 - 162*x^4*y - 243*x^4 + 6*x^2*y^5 - 27*x^2*y^3 + y^6 - 9*y^4 - 9*y^3;
p:=11;
N:=20;
data:=coleman_data(Q,p,N);
L,v:=effective_chabauty(data:bound:=1000,e:=75);

L;
/*
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(11^8),
        b := [ 1 + O(11^20), -3 + O(11^16), 9 + O(11^16), -27 + O(11^16), 27 + O(11^16), O(11^8) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(11^20),
        b := [ 1 + O(11^20), O(11^20), O(11^20), O(11^20), 27 + O(11^16), O(11^8) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(11^20),
        b := [ 1 + O(11^20), O(11^20), O(11^20), O(11^20), O(11^8), 11487432465893101 + O(11^16) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(11^9),
        b := [ 1 + O(11^20), 3 + O(11^20), 9 + O(11^20), 27 + O(11^20), -27 + O(11^9), -1178973805 + O(11^9) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(11^9),
        b := [ 1 + O(11^20), -3 + O(11^9), 9 + O(11^9), -27 + O(11^9), -81 + O(11^9), -1178973805 + O(11^9) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(11^18),
        b := [ 1 + O(11^20), -3 + O(11^18), 9 + O(11^18), -27 + O(11^18), -81 + O(11^9), 81 + O(11^9) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(11^9),
        b := [ 1 + O(11^20), 3 + O(11^20), 9 + O(11^20), 27 + O(11^20), 27 + O(11^9), -1178973805 + O(11^9) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(11^9),
        b := [ 1 + O(11^20), -3 + O(11^9), 9 + O(11^9), -27 + O(11^9), 81 + O(11^9), -1178973805 + O(11^9) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(11^18),
        b := [ 1 + O(11^20), -3 + O(11^18), 9 + O(11^18), -27 + O(11^18), 81 + O(11^9), 81 + O(11^9) ],
        inf := false>
]
*/

Q_points(data,1000);
/*
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(11^20),
        b := [ 1 + O(11^20), -3 + O(11^20), 9 + O(11^20), -27 + O(11^20), -81 + O(11^20), -336374997466280004560 + O(11^20) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(11^20),
        b := [ 1 + O(11^20), 3 + O(11^20), 9 + O(11^20), 27 + O(11^20), -27 + O(11^20), -336374997466280004560 + O(11^20) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(11^20),
        b := [ 1 + O(11^20), -3 + O(11^20), 9 + O(11^20), -27 + O(11^20), -81 + O(11^20), 81 + O(11^20) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(11^20),
        b := [ 1 + O(11^20), -3 + O(11^20), 9 + O(11^20), -27 + O(11^20), 81 + O(11^20), -336374997466280004560 + O(11^20) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(11^20),
        b := [ 1 + O(11^20), 3 + O(11^20), 9 + O(11^20), 27 + O(11^20), 27 + O(11^20), -336374997466280004560 + O(11^20) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(11^20),
        b := [ 1 + O(11^20), -3 + O(11^20), 9 + O(11^20), -27 + O(11^20), 81 + O(11^20), 81 + O(11^20) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(11^20),
        b := [ 1 + O(11^20), O(11^20), O(11^20), O(11^20), O(11^20), 168187498733140002361 + O(11^20) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(11^20),
        b := [ 1 + O(11^20), -3 + O(11^20), 9 + O(11^20), -27 + O(11^20), 27 + O(11^20), O(11^20) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(11^20),
        b := [ 1 + O(11^20), O(11^20), O(11^20), O(11^20), 27 + O(11^20), O(11^20) ],
        inf := true>
]
*/

v;
/*
[
    [ 1 + O(11^10), O(11^10), O(11^10), O(11^10), O(11^10), O(11^10), O(11^10) ],
    [ O(11^10), 1 + O(11^10), O(11^10), O(11^10), O(11^10), O(11^10), O(11^10) ],
    [ O(11^10), O(11^10), 1 + O(11^10), O(11^10), O(11^10), O(11^10), O(11^10) ],
    [ O(11^10), O(11^10), O(11^10), 1 + O(11^10), O(11^10), O(11^10), O(11^10) ],
    [ O(11^10), O(11^10), O(11^10), O(11^10), 1 + O(11^10), O(11^10), O(11^10) ],
    [ O(11^10), O(11^10), O(11^10), O(11^10), O(11^10), 1 + O(11^10), O(11^10) ],
    [ O(11^10), O(11^10), O(11^10), O(11^10), O(11^10), O(11^10), 1 + O(11^10) ]
]
*/

