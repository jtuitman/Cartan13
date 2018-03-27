SetPath("./Coleman-1.1");
load "coleman.m";
load "frobenius.m";
load "hodge.m";
load "heights.m";

Q:=y^4 + 5*x^4 - 6*x^2*y^2 + 6*x^3 + 26*x^2*y + 10*x*y^2 - 10*y^3 - 32*x^2 -40*x*y + 24*y^2 + 32*x - 16*y; // equation of the curve
p:=17;      // prime number p
N:=20;      // initial p-adic precision

Qp:=pAdicField(p,N);

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

FF:=fun_field(data);

bpt:=Zeros(FF.1)[1]; // point [0,0] as place on the function field


////////////////////
// list of points //
////////////////////

bound:=1000;
Qpoints:=Q_points(data,bound);
Qppoints:=Qp_points(data:points:=Qpoints); // first 2 points are infinite, last point finite bad, all other points good

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

Z1,A11:=hecke_corr(data,11,10:basis0:=basis0,basis1:=basis1);
eta1,betafil1,gammafil1:=hodge_data(data,denombasis,Z1,bpt); 

G1:=frob_struc(data,denombasis,Z1,eta1,[0,0]); // matrix of Frobenius structure on A_Z1(b)

G1_list:=[**]; // evaluations of G1 at Teichmuellers of all good points (0 if bad)
for i:=1 to #Qppoints do
  if is_bad(Qppoints[i],data) then
    G1_list[i]:=0;
  else
    P:=teichpoints[i];
    pt:=[IntegerRing()!(P`x),IntegerRing()!(P`b)[2]];
    G1_list[i]:=eval_mat_R(G1,pt,r);
  end if;
end for;

PhiAZ1b:=[**]; // Frobenius on the phi-modules A_Z1(b,P) (0 if P bad)
for i:=1 to #G1_list do
  if G1_list[i] ne 0 then
    PhiAZ1b[i]:=parallel_transport(teichpoints[i],Qppoints[i],denombasis,Z1,eta1,data:prec:=100)*frob_equiv_iso(G1_list[i],data);
  else
    PhiAZ1b[i]:=0;
  end if;
end for;

gammafil1_list:=[**]; // evaluations of gammafil1 at all good points (0 if bad)
for i:=1 to #G1_list do
  if G1_list[i] ne 0 then
    gammafil1_list[i]:=evalf0(ChangeRing(gammafil1,LaurentSeriesRing(BaseRing(gammafil1))),Qppoints[i],data);
  else
    gammafil1_list[i]:=0;
  end if;
end for;

////////////////////////////
// second correspondence: //
////////////////////////////

Z2,A7:=hecke_corr(data,7,10:basis0:=basis0,basis1:=basis1);
eta2,betafil2,gammafil2:=hodge_data(data,denombasis,Z2,bpt); 

G2:=frob_struc(data,denombasis,Z2,eta2,[0,0]); // matrix of Frobenius structure on A_Z2(b)
G2_list:=[**]; // evaluations of G2 at Teichmuellers of all good points (0 if bad)
for i:=1 to #Qppoints do
  if is_bad(Qppoints[i],data) then
    G2_list[i]:=0;
  else
    P:=teichpoints[i];
    pt:=[IntegerRing()!(P`x),IntegerRing()!(P`b)[2]];
    G2_list[i]:=eval_mat_R(G2,pt,r);
  end if;
end for;

PhiAZ2b:=[**]; // Frobenius on the phi-modules A_Z2(b,P) (0 if P bad)
for i:=1 to #G2_list do
  if G2_list[i] ne 0 then
    PhiAZ2b[i]:=parallel_transport(teichpoints[i],Qppoints[i],denombasis,Z2,eta2,data:prec:=100)*frob_equiv_iso(G2_list[i],data);
  else
    PhiAZ2b[i]:=0;
  end if;
end for;

gammafil2_list:=[**]; // evaluations of gammafil2 at Teichmuellers of all good points (0 if bad)
for i:=1 to #G2_list do
  if G2_list[i] ne 0 then
    gammafil2_list[i]:=evalf0(ChangeRing(gammafil2,LaurentSeriesRing(BaseRing(gammafil2))),Qppoints[i],data);
  else
    gammafil2_list[i]:=0;
  end if;
end for;
