load "frobenius.m";

Q:=y^4 + 5*x^4 - 6*x^2*y^2 + 6*x^3 + 26*x^2*y + 10*x*y^2 - 10*y^3 - 32*x^2 -40*x*y + 24*y^2 + 32*x - 16*y; // equation of the curve
p:=17;    // prime number p
N:=20;    // initial p-adic precision
b:=[0,0]; // base point

r,Delta,s:=auxpolys(Q);

// put in a basis omega[i]dx/(dQ/dy) for H^1(X) by hand:

omega:=[Zxy|];
omega[1]:=3;
omega[2]:=3*x;
omega[3]:=3*y;
omega[4]:=-160*x^4+736*x^3-16*x^2*y+436*x^2-440*x*y+68*y^2;
omega[5]:=-80*x^3+132*x^2-40*x*y+68*y^2-96;
omega[6]:=-48*x^2*y+84*x^2+216*x*y-12*y^2-160*x+272;

denom:=3; // the omega's have to be integral, but a denominator can be specified to make the basis symplectic

///////////////////////////
// first correspondence: //
///////////////////////////

// Z1:=hecke_corr(Q,11,10,omega,denom);
Z1:=Matrix(RationalField(),6,6,[[0,-976,-1104,10,-6,18],[976,0,-816,-3,1,3],[1104,816,0,-3,3,-11],[-10,3,3,0,0,0],[6,-1,-3,0,0,0],[-18,-3,11,0,0,0]]);
eta1:=-(132*x^2+148*x*y+24*y^2);

Phi1:=frob_struc(Q,p,N,Z1,eta1,b,omega,denom);
pt:=[0,0];
Phi1_at_pt:=eval_mat_R(Phi1,pt,r);

////////////////////////////
// second correspondence: //
////////////////////////////

// Z2:=hecke_corr(Q,7,10,omega,denom);
Z2:=Matrix(RationalField(),6,6,[[0,112,-656,-6,6,6],[-112,0,-2576,15,9,27],[656,2576,0,3,3,-3],[6,-15,-3,0,0,0],[-6,-9,-3,0,0,0],[-6,-27,3,0,0,0]]);
eta2:=3*(-40*x^2+148*x*y+36*y^2);

Phi2:=frob_struc(Q,p,N,Z2,eta2,b,omega,denom);
pt:=[0,0];
Phi2_at_pt:=eval_mat_R(Phi2,pt,r);
