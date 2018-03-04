load "frobenius.m";
load "hodge.m";

Q:=y^4 + 5*x^4 - 6*x^2*y^2 + 6*x^3 + 26*x^2*y + 10*x*y^2 - 10*y^3 - 32*x^2 -40*x*y + 24*y^2 + 32*x - 16*y; // equation of the curve
p:=17;      // prime number p
N:=20;      // initial p-adic precision
bpt:=[0,0]; // base point

r,Delta,s:=auxpolys(Q);

// put in a basis omega[i]dx/(dQ/dy) for H^1(Y) by hand:

omega:=[Zxy|];
omega[1]:=3;
omega[2]:=3*x;
omega[3]:=3*y;
omega[4]:=-160*x^4+736*x^3-16*x^2*y+436*x^2-440*x*y+68*y^2;
omega[5]:=-80*x^3+132*x^2-40*x*y+68*y^2-96;
omega[6]:=-48*x^2*y+84*x^2+216*x*y-12*y^2-160*x+272;
omega[7]:=3*x^2;
omega[8]:=3*x*y;
omega[9]:=3*y^2;

denomomega:=3; // the omega's have to be integral, but a denominator can be specified

basis0:=[]; // first kind
for i:=1 to 3 do
  basis0[i]:=Coefficients(reduce_mod_Q_exact(omega[i]*s,Q));
end for;

basis1:=[]; // second kind
for i:=1 to 3 do
  basis1[i]:=Coefficients(reduce_mod_Q_exact(omega[i+3]*s,Q));
end for;

basis2:=[]; // basis for H^1(Y) over H^1(X)
for i:=1 to 3 do
  basis2[i]:=Coefficients(reduce_mod_Q_exact(omega[i+6]*s,Q));
end for;

denombasis:=denomomega*LeadingCoefficient(Delta);

// b0 cat b1 is the basis for H^1(X) given by omega[i]*dx/(dQ/dy), multiplied by denom*LeadingCoefficient(Delta).

data:=coleman_data(Q,p,N:useU:=true,basis0:=basis0,basis1:=basis1,basis2:=basis2);

////////////////////
// list of points //
////////////////////

bound:=1000;
Qpoints:=Q_points(data,bound);

Qppoints:=Qp_points(data); // first 2 points are infinite, last point finite bad, all other points good

teichpoints:=[**]; // compute Teichmueller representatives of good points
for i:=1 to #Qppoints do
  if is_bad(Qppoints[i],data) then
    teichpoints[i]:=0;
  else
    teichpoints[i]:=teichmueller_pt(Qppoints[i],data);
  end if;
end for;

///////////////////////////
// first correspondence: //
///////////////////////////

// Z1:=hecke_corr(data,11,10:basis0:=basis0,basis1:=basis1);
Z1:=Matrix(RationalField(),6,6,[[0,-976,-1104,10,-6,18],[976,0,-816,-3,1,3],[1104,816,0,-3,3,-11],[-10,3,3,0,0,0],[6,-1,-3,0,0,0],[-18,-3,11,0,0,0]]);
eta1:=-(132*x^2+148*x*y+24*y^2);
// hodge_data(data,denombasis,Z1); 

G1:=frob_struc(data,Z1,eta1,bpt,denombasis);
G1_list:=[**];
for i:=1 to #Qppoints do
  if is_bad(Qppoints[i],data) then
    G1_list[i]:=0;
  else
    P:=teichpoints[i];
    pt:=[IntegerRing()!(P`x),IntegerRing()!(P`b)[2]];
    G1_list[i]:=eval_mat_R(G1,pt,r);
  end if;
end for;

////////////////////////////
// second correspondence: //
////////////////////////////

// Z2:=hecke_corr(data,7,10:basis0:=basis0,basis1:=basis1);
Z2:=Matrix(RationalField(),6,6,[[0,112,-656,-6,6,6],[-112,0,-2576,15,9,27],[656,2576,0,3,3,-3],[6,-15,-3,0,0,0],[-6,-9,-3,0,0,0],[-6,-27,3,0,0,0]]);
eta2:=3*(-40*x^2+148*x*y+36*y^2);
// hodge_data(data,denombasis,Z2); 

G2:=frob_struc(data,Z2,eta2,bpt,denombasis);
G2_list:=[**];
for i:=1 to #Qppoints do
  if is_bad(Qppoints[i],data) then
    G2_list[i]:=0;
  else
    P:=teichpoints[i];
    pt:=[IntegerRing()!(P`x),IntegerRing()!(P`b)[2]];
    G2_list[i]:=eval_mat_R(G2,pt,r);
  end if;
end for;
