SetPath("./Coleman-1.3");

load "coleman.m";
load "frobenius.m";
load "hodge.m";
load "heights.m";

Q:=-16*x^3*y+32*x^3+24*x^2*y^2-40*x^2*y-32*x^2-10*x*y^3+10*x*y^2+26*x*y+6*x+y^4-6*y^2+5; // x^4*Evaluate(Q,[1/x,y/x])
p:=17;      // prime number p
N:=20;      // initial p-adic precision
prec:=25;   // t-adic precision used in expansions

Qp:=pAdicField(p,N);

r,Delta,s:=auxpolys(Q);
lc:=LeadingCoefficient(Delta);

// put in a basis omega[i]dx/z for H^1(Y) by hand:

omega:=[Qxy|];
omega[1]:=-x/lc;
omega[2]:=-1/lc;
omega[3]:=-y/lc;
omega[4]:=(768/5*x^2*y-448/5*x*y^2-1536/5*x^2+96*x*y+16*y^2+2272/15*x-1648/15*y+1712/15)/lc;
omega[5]:=(128/7*x^2*y^2-5056/35*x^2*y+576/35*x*y^2+7552/35*x^2-816/7*x*y+136/7*y^2+10736/105*x-1072/15*y-184/105)/lc;
omega[6]:=(-448/5*x^2*y+288/5*x*y^2+896/5*x^2-80*x*y-8*y^2-2272/15*x+96/5*y-1432/15)/lc; 
omega[7]:=x^2/lc;
omega[8]:=x*y/lc;
omega[9]:=y^2/lc;

// p should not be 3,5,7 because of denominator in basis

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

// basis0 cat basis1 is the basis for H^1(X) given by omega[i]*dx/z.

data:=coleman_data(Q,p,N:useU:=true,basis0:=basis0,basis1:=basis1,basis2:=basis2);

FF:=fun_field(data);

bpt:=Zeros(FF.1+1)[1]; // point [0,-1] as place on the function field


////////////////////
// list of points //
////////////////////

bound:=1000;
Qpoints:=Q_points(data,bound);
Qppoints:=Qp_points(data:points:=Qpoints); // first 4 points are infinite, last point finite bad, all other points good
numberofpoints:=#Qppoints;

teichpoints:=[**]; // compute Teichmueller representatives of good points
for i:=1 to numberofpoints do
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
eta1,betafil1,gammafil1:=hodge_data(data,Z1,bpt); 

G1:=frob_struc(data,Z1,eta1,[0,-1]); // matrix of Frobenius structure on A_Z1(b)
G1_list:=[**]; // evaluations of G1 at Teichmuellers of all good points (0 if bad)
for i:=1 to numberofpoints do
  if is_bad(Qppoints[i],data) then
    G1_list[i]:=0;
  else
    P:=teichpoints[i];
    pt:=[IntegerRing()!(P`x),IntegerRing()!(P`b)[2]];
    G1_list[i]:=eval_mat_R(G1,pt,r);
  end if;
end for;

PhiAZ1b:=[**]; // Frobenius on the phi-modules A_Z1(b,P) (0 if P bad)
for i:=1 to numberofpoints do
  if G1_list[i] ne 0 then
    PhiAZ1b[i]:=parallel_transport(teichpoints[i],Qppoints[i],Z1,eta1,data:prec:=prec)*frob_equiv_iso(G1_list[i],data);
  else
    PhiAZ1b[i]:=0;
  end if;
end for;

PhiAZ1b_to_z:=[**]; // Frobenius on the phi-modules A_Z1(b,z) for z in residue disk of P (0 if P bad)
for i:=1 to numberofpoints do
  if G1_list[i] ne 0 then
    PhiAZ1b_to_z[i]:=parallel_transport_to_z(Qppoints[i],Z1,eta1,data:prec:=prec)*PhiAZ1b[i];
  else
    PhiAZ1b_to_z[i]:=0;
  end if;
end for;

gammafil1_list:=[**]; // evaluations of gammafil1 at all good points (0 if bad)
gammafil1_listb_to_z:=[**]; // evaluations of gammafil1 at all good points (0 if bad)
for i:=1 to numberofpoints do
    gammafil1_listb_to_z[i]:=evalf0_bad(Qppoints[i],gammafil1,data,N,prec);
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
eta2,betafil2,gammafil2:=hodge_data(data,Z2,bpt); 

G2:=frob_struc(data,Z2,eta2,[0,-1]); // matrix of Frobenius structure on A_Z2(b)
G2_list:=[**]; // evaluations of G2 at Teichmuellers of all good points (0 if bad)
for i:=1 to numberofpoints do
  if is_bad(Qppoints[i],data) then
    G2_list[i]:=0;
  else
    P:=teichpoints[i];
    pt:=[IntegerRing()!(P`x),IntegerRing()!(P`b)[2]];
    G2_list[i]:=eval_mat_R(G2,pt,r);
  end if;
end for;

PhiAZ2b:=[**]; // Frobenius on the phi-modules A_Z2(b,P) (0 if P bad)
for i:=1 to numberofpoints do
  if G2_list[i] ne 0 then
    PhiAZ2b[i]:=parallel_transport(teichpoints[i],Qppoints[i],Z2,eta2,data:prec:=prec)*frob_equiv_iso(G2_list[i],data);
  else
    PhiAZ2b[i]:=0;
  end if;
end for;

PhiAZ2b_to_z:=[**]; // Frobenius on the phi-modules A_Z2(b,z) for z in residue disk of P (0 if P bad)
for i:=1 to numberofpoints do
  if G2_list[i] ne 0 then
    PhiAZ2b_to_z[i]:=parallel_transport_to_z(Qppoints[i],Z2,eta2,data:prec:=prec)*PhiAZ2b[i];
  else
    PhiAZ2b_to_z[i]:=0;
  end if;
end for;

gammafil2_list:=[**]; // evaluations of gammafil2 at Teichmuellers of all good points (0 if bad)
gammafil2_listb_to_z:=[**];
for i:=1 to numberofpoints do
  gammafil2_listb_to_z[i]:=evalf0_bad(Qppoints[i],gammafil2,data,N,prec);
  if G2_list[i] ne 0 then
    gammafil2_list[i]:=evalf0(ChangeRing(gammafil2,LaurentSeriesRing(BaseRing(gammafil2))),Qppoints[i],data);
  else
    gammafil2_list[i]:=0;
  end if;
end for;


/////////////
// heights //
/////////////

P1:=Qppoints[19];
P3:=Qppoints[9];
P4:=Qppoints[6];
P6:=Qppoints[5]; // base point

eqsplit:=Matrix(RationalField(),6,3,[ 1, 0, 0, 0, 1, 0, 0, 0, 1, 224/3, -880/3, 0, -880/3, -1696/3, 0, 0, 0, 0 ]); // equivariant splitting of Hodge filtration, put in by hand for now
height1_P1:=height(PhiAZ1b[19],betafil1,gammafil1_list[19],eqsplit,data); 
height1_P3:=height(PhiAZ1b[9],betafil1,gammafil1_list[9],eqsplit,data); 
height1_P4:=height(PhiAZ1b[6],betafil1,gammafil1_list[6],eqsplit,data);

q:=3;

_,Aq:=hecke_corr(data,q,20:basis0:=basis0,basis1:=basis1);                   // Hecke operator at q on H^1_dR
Aq_small:=ExtractBlock(Aq,1,1,3,3);                                          // Hecke operator at q on H^0(Omega^1), A3 is wrong because of denominator 3 in basis, but A3_small is not affected
m:=CharacteristicPolynomial(Aq_small);

assert IsZero(eqsplit*ExtractBlock(Aq,4,4,3,3) - Transpose(Aq)*eqsplit);     // Test equivariant splitting of Hodge filtration

E1P1:=Vector(Qp,3,[PhiAZ1b[19][i+1,1]:i in [1..3]]);                         // AJ_b(P1) which generates H^0(Omega^1)^* over K                               
basisH0star:=[];
for i:=0 to 2 do
  basisH0star:=Append(basisH0star,Eltseq(E1P1*(ChangeRing(Aq_small,Qp)^i))); // basis for H^0(Omega^1)^* generated by powers of iota(aq) acting on E1P1
end for; 

E1_E2_P1 := E1_tensor_E2(PhiAZ1b[19],betafil1,basisH0star,m,data); // P1
E1_E2_P3 := E1_tensor_E2(PhiAZ1b[9],betafil1,basisH0star,m,data);  // P3
E1_E2_P4 := E1_tensor_E2(PhiAZ1b[6],betafil1,basisH0star,m,data);  // P4

Nend:=Floor(N/2);
Qp:=pAdicField(p,Nend); // TODO For now assuming floor(N/2) digits are correct
S:=LaurentSeriesRing(Qp,prec);

F1_list:=[**];
for i:=1 to numberofpoints do
  if G1_list[i] eq 0 then
    F1_list[i]:=0;
  else
    T1:=ZeroMatrix(S,4,4);
    T1[1,1]:=height(PhiAZ1b_to_z[i],betafil1,gammafil1_listb_to_z[i],eqsplit,data);
    for j:=2 to 4 do
      T1[1,j]:=Eltseq(E1_tensor_E2(PhiAZ1b_to_z[i],betafil1,basisH0star,m,data))[j-1];
    end for;
    T1[2,1]:=height1_P1;
    T1[3,1]:=height1_P3;
    T1[4,1]:=height1_P4; 
    for j:=2 to 4 do
      T1[2,j]:=Eltseq(E1_E2_P1)[j-1];
      T1[3,j]:=Eltseq(E1_E2_P3)[j-1];
      T1[4,j]:=Eltseq(E1_E2_P4)[j-1];
    end for;
    F1_list[i]:=Determinant(T1);
  end if;
end for;

F2_list:=[**];
for i:=1 to numberofpoints do
  if G2_list[i] eq 0 then
    F2_list[i]:=0;
  else
    T2:=ZeroMatrix(S,4,4);
    T2[1,1]:=height(PhiAZ2b_to_z[i],betafil2,gammafil2_listb_to_z[i],eqsplit,data);
    for j:=2 to 4 do
      T2[1,j]:=Eltseq(E1_tensor_E2(PhiAZ2b_to_z[i],betafil2,basisH0star,m,data))[j-1];
    end for;
    T2[2,1]:=height1_P1;
    T2[3,1]:=height1_P3;
    T2[4,1]:=height1_P4; 
    for j:=2 to 4 do
      T2[2,j]:=Eltseq(E1_E2_P1)[j-1];
      T2[3,j]:=Eltseq(E1_E2_P3)[j-1];
      T2[4,j]:=Eltseq(E1_E2_P4)[j-1];
    end for;
    F2_list[i]:=Determinant(T2);
  end if;
end for;        


///////////////////////
// test known points //
///////////////////////

for i in [5,6,9,19] do
  assert Valuation(Evaluate(F1_list[i],0)) ge Nend;
  assert Valuation(Evaluate(F2_list[i],0)) ge Nend;
end for;


////////////////
// find zeros //
////////////////

Qptt:=PowerSeriesRing(Qp);
Zp:=pAdicRing(p,Precision(Qp));
Zpt:=PolynomialRing(Zp);

zero1_list:=[**];
for i:=1 to numberofpoints do
  if G1_list[i] eq 0 then
    zero1_list[i]:=0;
  else
    f:=F1_list[i];
    f:=Evaluate(Qptt!f,p*Qptt.1);
    val:=Valuation(f);
    f:=Zpt.1^val*(Zpt![Zp!c : c in Coefficients(f)]);
    zero1_list[i]:=my_roots_Zpt(f);
  end if;
end for;

zero2_list:=[**];
for i:=1 to numberofpoints do
  if G2_list[i] eq 0 then
    zero2_list[i]:=0;
  else
    f:=F2_list[i];
    f:=Evaluate(Qptt!f,p*Qptt.1);
    val:=Valuation(f);
    f:=Zpt.1^val*(Zpt![Zp!c : c in Coefficients(f)]);
    zero2_list[i]:=my_roots_Zpt(f);
  end if;
end for;
