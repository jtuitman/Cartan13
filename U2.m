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
omega[6]:=(-448/5*x^2*y+288/5*x*y^2+896/5*x^2-80*x*y-8*y^2-2272/15*x+96/5*y-1432/5)/lc;
omega[7]:=x^2/lc;
omega[8]:=x*y/lc;
omega[9]:=y^2/lc;

// p should not be 3,5,7 because of this choice of basis

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

Z1,A11:=hecke_corr(data,11,10:basis0:=basis0,basis1:=basis1);
eta1,betafil1,gammafil1:=hodge_data(data,Z1,bpt);

Z2,A7:=hecke_corr(data,7,10:basis0:=basis0,basis1:=basis1);
eta2,betafil2,gammafil2:=hodge_data(data,Z2,bpt);

print "eta1,betafil1,gammafil1","/n",eta1,betafil1,gammafil1;
print "eta2,betafil2,gammafil2","/n",eta2,betafil2,gammafil2;

