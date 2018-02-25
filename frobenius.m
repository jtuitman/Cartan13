//////////////////////////////////////////////////
// Functions for computing Frobenius structures //
//////////////////////////////////////////////////


SetPath("./Coleman-1.1");
load "coleman.m";


hecke_corr:=function(data,q,N:b0:=[],b1:=[])

  // compute the matrix of the correspondence Z_q constructed from 
  // the Hecke operator A_q w.r.t. the basis of H^1(X) given by [b0,b1].
  // assumes that this basis is symplectic

  Q:=data`Q; g:=data`g;

  data:=coleman_data(Q,q,N:b0:=b0,b1:=b1);

  F:=data`F;
  Aq:=Transpose(F)+q*Transpose(F)^(-1);  
  C:=ZeroMatrix(RationalField(),2*g,2*g);
  for i:=1 to g do
    C[i,g+i]:=-1;
  end for;
  for i:=1 to g do 
    C[g+i,i]:=1; 
  end for;
  Z:=(2*g*Aq-Trace(Aq)*IdentityMatrix(RationalField(),2*g))*C^(-1);
  
  for i:=1 to 2*g do
    for j:=1 to 2*g do
      Z[i,j]:=reduce_mod_pN_Q(Z[i,j],q,Floor(N/2)); // assumes that Z mod q^(N/2) will recover Z exactly, not rigorous for now
    end for;
  end for;

  return Z;

end function;


frob_struc:=function(data,Z,eta,bpt,omega,denom)

  // compute the matrix of the Frobenius structure Phi on A_Z using the basis for 
  // H^1(X) given by [omega, denom], the matrix Z of a nice correspondence and a
  // 1-form eta w.r.t. this basis.

  Q:=data`Q; p:=data`p; N:=data`N; g:=data`g; W0:=data`W0; Winf:=data`Winf; f0list:=data`f0list; finflist:=data`finflist; fendlist:=data`fendlist; 
  FH1U:=data`F; Nmax:=data`Nmax; basis:=data`basis; G0:=data`G0; Ginf:=data`Ginf; red_list_fin:=data`red_list_fin; red_list_inf:=data`red_list_inf; 
  basis:=data`basis; integrals:=data`integrals; quo_map:=data`quo_map; r:=data`r; Delta:=data`Delta; s:=data`s;

  d:=Degree(Q); lc:=LeadingCoefficient(Delta);

  O,Ox,S,R:=getrings(p,data`Nmax); // O = IntegerRing(p^Nmax), Ox = O[x], S = Ox[z,1/z], R = S[y]

  // Coerce Q into R:

  C:=[];
  for i:=1 to Degree(Q)+1 do
    C[i]:=Ox!(Coefficient(Q,i-1));
  end for;
  QR:=(R!C);

  // Coerce z = r/LeadingCoefficient(r) into R:

  zR:=(Ox!r)/LeadingCoefficient(r);

  // Coerce the f_i = f_{i,0}+f_{i,inf}+f_{i,end} into R:

  fRlist:=[];
  for i:=1 to 2*g do
    fRlist[i]:=Qxzzinvd_to_R(compute_F(Q,W0,Winf,f0list[i],finflist[i],fendlist[i]),Q,p,r,R,W0);
    fRlist[i]:=reduce_mod_Q(fRlist[i],QR,zR);
    fRlist[i]:=fRlist[i]-eval_R(fRlist[i],bpt,r); // make sure that f_i(bpt) = 0 
  end for;

  // The matrix of Frobenius on H^1(X) is the 6x6 top left corner of the matrix of Frobenius on H^1(U):

  FH1X:=ZeroMatrix(RationalField(),6,6);
  for i:=1 to 2*g do
    for j:=1 to 2*g do
      FH1X[i,j]:=FH1U[i,j];
    end for;
  end for;

  // Compute g0:

  A:=-Transpose(FH1X)*Z;
  A:=ChangeRing(A,R);
  g0:=[];
  for i:=1 to 2*g do
    g0[i]:=R!0;
    for j:=1 to 2*g do
      g0[i]:=g0[i]+A[i,j]*fRlist[j]; // multiplied by denom*lc.
    end for;
  end for;

  // Coerce s into R:

  C:=[];
  for i:=1 to Degree(s)+1 do
    C[i]:=Ox!(Coefficient(s,i-1));
  end for;
  sR:=(R!C);

  // Coerce basis elements of H^1(X) into R:

  basisR:=[];
  for i:=1 to 2*g do
    basisR[i]:=reduce_mod_Q((R![S.1^0*(Ox!c) : c in Coefficients(omega[i])])*sR,QR,zR);  
  end for;

  // Compute g0*omega:

  g0omega:=R!0;
  for i:=1 to 2*g do
    g0omega:=g0omega+g0[i]*basisR[i];
  end for;
  g0omega:=reduce_mod_Q(g0omega,QR,zR);
  g0omega:=Vector(Coefficients(g0omega));

  // Precompute Fp(y^i/r) for i=0,..,3 (in fact, this has already happened inside coleman_data, but is not available here)

  frobmatb0r:=froblift(Q,p,Nmax-1,r,Delta,s,W0);

  // Compute phi^(*)(eta)-p(eta):

  eta:=eta*s; // multiplied by denom*lc
  eta:=reduce_mod_Q_exact(eta,Q);
  eta:=Vector(Coefficients(eta));
  phi_eta:=frobenius(eta,Q,p,Nmax,r,frobmatb0r); // phi^(*)(eta)
  for i:=1 to d do
    phi_eta[i]:=phi_eta[i]-p*(S!Evaluate(eta[i],Ox.1)); // phi^(*)(eta)-p(eta), as vector in S^4
  end for;

  // Compute phi^*(omega):

  phiomega:=[];
  for i:=1 to 2*g do
    phiomegai:=frobenius(basis[i],Q,p,Nmax,r,frobmatb0r);
    phiomega[i]:=R!0;
    for j:=1 to d do
      phiomega[i]:=phiomega[i]+phiomegai[j]*R.1^(j-1);
    end for;
  end for;

  // Compute Z*f:

  Zf:=[];
  for i:=1 to 2*g do
    Zf[i]:=R!0;
    for j:=1 to 2*g do
      Zf[i]:=Zf[i]+ChangeRing(Z,R)[i,j]*fRlist[j]; // multiplied by denom*lc 
    end for;
  end for;

  // Compute phi*omega*Zf:

  phiomegaZf:=R!0;
  for i:=1 to 2*g do
    phiomegaZf:=reduce_mod_Q(phiomegaZf+phiomega[i]*Zf[i],QR,zR);
  end for;
  phiomegaZf:=Vector(Coefficients(phiomegaZf));

  // Compute c and h:

  sum:=phiomegaZf+denom*lc*phi_eta-g0omega; // multiplied by (denom*lc)^2 
  sum:=convert_to_Qxzzinvd(sum,Q);
  c,h0,hinf,hend:=reduce_with_fs(sum,Q,p,N,Nmax,r,W0,Winf,G0,Ginf,red_list_fin,red_list_inf,basis,integrals,quo_map); // multiplied by (denom*lc)^2 
  hR:=Qxzzinvd_to_R(compute_F(Q,W0,Winf,h0,hinf,hend),Q,p,r,R,W0);
  hR:=reduce_mod_Q(hR,QR,zR);
  hR:=hR-eval_R(hR,bpt,r); // make sure that h(bpt) = 0 

  // Correct for the factors denom*lc:

  f:=[];
  for i:=1 to 2*g do
    f[i]:=fRlist[i]/(denom*lc); // renormalised
  end for;

  g1:=[];
  for i:=1 to 2*g do
    g1[i]:=g0[i]/(denom*lc)+(O!c[i])/(denom*lc); // renormalised
  end for;

  h:=hR/(denom*lc)^2; // renormalised

  // Compute matrix G:

  G:=ZeroMatrix(R,2*g+2,2*g+2);
  G[1,1]:=1;
  for i:=1 to 2*g do
    G[i+1,1]:=f[i];
  end for;
  G[8,1]:=h;
  for i:=1 to 2*g do
    for j:=1 to 2*g do
      G[i+1,j+1]:=FH1X[i,j];
    end for;
  end for;
  for i:=1 to 2*g do
    G[2*g+2,i+1]:=g1[i];
  end for;
  G[2*g+2,2*g+2]:=p;

  return G;
end function;

