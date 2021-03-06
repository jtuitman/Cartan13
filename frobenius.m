//////////////////////////////////////////////////
// Functions for computing Frobenius structures //
//////////////////////////////////////////////////


hecke_corr:=function(data,q,N:basis0:=[],basis1:=[])

  // compute the matrix of the correspondence Z_q constructed from 
  // the Hecke operator A_q w.r.t. the basis of H^1(X) given by [b0,b1].
  // assumes that this basis is symplectic

  Q:=data`Q; g:=data`g; d:=Degree(Q);

  if basis0 ne [] then
    v0:=Minimum(&cat[[Valuation(coef,q):coef in &cat[Coefficients(basis0[i][j]):j in [1..d]]]: i in [1..g]]); // valuation basis0
  else
    v0:=0;
  end if;

  if basis1 ne [] then
    v1:=Minimum(&cat[[Valuation(coef,q):coef in &cat[Coefficients(basis1[i][j]):j in [1..d]]]: i in [1..g]]); // valuation basis1
  else
    v1:=0;
  end if;

  v:=Minimum([v0,v1]);

  // multiply by constant to remove denominator basis0 and basis1

  if v lt 0 then
    for i:=1 to g do
      for j:=1 to d do
        basis0[i][j]:=q^(-v)*basis0[i][j];
        basis1[i][j]:=q^(-v)*basis1[i][j];
      end for;
    end for;
  end if;

  data:=coleman_data(Q,q,N:basis0:=basis0,basis1:=basis1);

  F:=data`F;
  Aq:=Transpose(F)+q*Transpose(F)^(-1);  
  Aqtimes2g:=2*g*Aq;
  C:=ZeroMatrix(RationalField(),2*g,2*g);
  for i:=1 to g do
    C[i,g+i]:=-1;
  end for;
  for i:=1 to g do 
    C[g+i,i]:=1; 
  end for;
  Z:=(Aqtimes2g-Trace(Aq)*IdentityMatrix(RationalField(),2*g))*C^(-1);

  // assume that precision q^(N/2) is sufficient to recover matrices Z and 2g*Aq exactly, not rigorous for now
  
  for i:=1 to 2*g do
    for j:=1 to 2*g do
      Aqtimes2g[i,j]:=reduce_mod_pN_Q(Aqtimes2g[i,j],q,Floor(N/2)); 
      Z[i,j]:=reduce_mod_pN_Q(Z[i,j],q,Floor(N/2));      
    end for;
  end for;

  Aq:=Aqtimes2g/(2*g);

  return Z, Aq;

end function;


frob_struc:=function(data,Z,eta,bpt)

  // Compute the matrix G of the (inverse) Frobenius structure on A_Z.

  Q:=data`Q; p:=data`p; N:=data`N; g:=data`g; W0:=data`W0; Winf:=data`Winf; f0list:=data`f0list; finflist:=data`finflist; fendlist:=data`fendlist; 
  FH1U:=data`F; Nmax:=data`Nmax; basis:=data`basis; G0:=data`G0; Ginf:=data`Ginf; red_list_fin:=data`red_list_fin; red_list_inf:=data`red_list_inf; 
  basis:=data`basis; integrals:=data`integrals; quo_map:=data`quo_map; r:=data`r; Delta:=data`Delta; s:=data`s;

  d:=Degree(Q); lc:=LeadingCoefficient(Delta);

  O,Ox,S,R:=getrings(p,Nmax); // O = IntegerRing(p^Nmax), Ox = O[x], S = Ox[z,1/z], R = S[y]

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

  // The matrix of Frobenius on H^1(X) is the 2gx2g top left corner of the matrix of Frobenius on H^1(U):

  FH1X:=ZeroMatrix(RationalField(),2*g,2*g);
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
      g0[i]:=g0[i]+A[i,j]*fRlist[j]; 
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
    basisR[i]:=R![S.1^0*(Ox!c) :c in Eltseq(basis[i])];
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

  // determine the number of points at infinity

  FF:=fun_field(data);
  infplaces:=InfinitePlaces(FF);
  infpoints:=0;
  for i:=1 to #infplaces do
    infpoints:=infpoints+Degree(infplaces[i]);
  end for;

  // Compute phi^(*)(eta)-p(eta):

  eta_new:=[];
  for i:=1 to d do
    sum:=0;
    for j:=1 to infpoints-1 do
      sum:=sum+Qx!(eta[j]*(PolynomialRing(RationalField())!basis[2*g+j][i])); // fix this
    end for;
    eta_new[i]:=sum;
  end for;
  eta:=eta_new;

  phi_eta:=frobenius(eta,Q,p,Nmax,r,frobmatb0r); // phi^(*)(eta)
  for i:=1 to d do
    phi_eta[i]:=phi_eta[i]-p*(S!Ox!eta[i]); // phi^(*)(eta)-p(eta), as vector in S^4
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
      Zf[i]:=Zf[i]+ChangeRing(Z,R)[i,j]*fRlist[j]; 
    end for;
  end for;

  // Compute phi*omega*Zf:

  phiomegaZf:=R!0;
  for i:=1 to 2*g do
    phiomegaZf:=reduce_mod_Q(phiomegaZf+phiomega[i]*Zf[i],QR,zR);
  end for;
  phiomegaZf:=Vector(Coefficients(phiomegaZf));

  // Compute c and h:

  sum:=phiomegaZf+phi_eta-g0omega;
  sum:=convert_to_Qxzzinvd(sum,Q);
  c,h0,hinf,hend:=reduce_with_fs(sum,Q,p,N,Nmax,r,W0,Winf,G0,Ginf,red_list_fin,red_list_inf,basis,integrals,quo_map); 
  hR:=Qxzzinvd_to_R(compute_F(Q,W0,Winf,h0,hinf,hend),Q,p,r,R,W0);
  hR:=reduce_mod_Q(hR,QR,zR);
  hR:=hR-eval_R(hR,bpt,r); // make sure that h(bpt) = 0 

  g1:=[];
  for i:=1 to 2*g do
    g1[i]:=g0[i]+(O!c[i]);
  end for;

  // Compute matrix G:

  G:=ZeroMatrix(R,2*g+2,2*g+2);
  G[1,1]:=1;
  for i:=1 to 2*g do
    G[i+1,1]:=fRlist[i];
  end for;
  G[2*g+2,1]:=hR;
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

