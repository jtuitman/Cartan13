/////////////////////////////////////
// functions for computing heights //
/////////////////////////////////////


function frob_equiv_iso(G,data)

  // Given a matrix G defining Frob^-1 on a mixed extension this returns s^phi.

  p:=data`p; N:=data`N; g:=data`g; F:=data`F;

  Qp:=pAdicField(p,N);  

  FQp:=ZeroMatrix(Qp,2*g,2*g);
  for i:=1 to 2*g do
    for j:=1 to 2*g do
      FQp[i,j]:=Qp!F[i,j];
    end for;
  end for;
  
  fx:=Matrix(Qp,2*g,1,[Qp!G[i+1,1]:i in [1..2*g]]);
  gx:=Matrix(Qp,1,2*g,[Qp!G[2*g+2,i+1]:i in [1..2*g]]);
  hx:=Matrix(Qp,1,1,[Qp!G[2*g+2,1]]);
  
  I := IdentityMatrix(Qp,2*g);
  
  hxprime := 1/(1-p)*(gx*((I-FQp)^-1)*fx+hx);
  
  s_phi := IdentityMatrix(Qp,2*g+2);
  for i in [1..2*g] do
    s_phi[i+1,1] := ((I-FQp)^(-1)*fx)[i,1];
    s_phi[2*g+2,i+1] := (gx*(FQp-p*I)^-1)[1,i];
  end for;
  s_phi[2*g+2,1] := hxprime[1,1];

  return s_phi;

end function;


function parallel_transport(P1,P2,denombasis,Z,eta,data:prec:=0)

  // Computes the parallel transport map of the unipotent 
  // connection Lambda defined by Z,eta from P1 to P2.

  // TODO: on a bad disk, make sure to take local 
  // parameters centered at the (very) bad point.

  x1:=P1`x; x2:=P2`x; b1:=P1`b; b2:=P2`b; Q:=data`Q; p:=data`p; N:=data`N; W0:=data`W0; Winf:=data`Winf; r:=data`r; basis:=data`basis; g:=data`g;
  d:=Degree(Q); lc_r:=LeadingCoefficient(r); K:=Parent(x1);

  if not lie_in_same_disk(P1,P2,data) then
    error "the points do not lie in the same residue disk";
  end if;

  xt,bt,index:=local_coord(P1,prec,data); 

  if index eq 0 then // x-x(P1) is the local coordinate
    val:=Valuation(x2-x1);
  else // b[index]-b[index](P1) is the local coordinate
    val:=Valuation(b2[index]-b1[index]);
  end if;

  Qt<t>:=LaurentSeriesRing(RationalField(),prec);
  xt:=Qt!xt;
  btnew:=[Qt|];
  for i:=1 to d do
    btnew[i]:=Qt!bt[i];
  end for;
  bt:=Vector(btnew);

  Kt:=LaurentSeriesRing(K,prec);
  OK:=RingOfIntegers(K);
  OKt:=LaurentSeriesRing(OK,prec);

  if P1`inf then
    xt:=1/xt;
    xt:=Qt!Kt!xt; 
    Winv:=W0*Winf^(-1);          
    bt:=bt*Transpose(Evaluate(Winv,xt));
    for i:=1 to d do
      bt[i]:=Qt!(Kt!bt[i]);
    end for; 
  end if;

  if P1`inf or not is_bad(P1,data) then 
    denom:=Qt!Kt!(1/Evaluate(r,xt));
  else
    Qp:=pAdicField(p,N);
    Qpx:=PolynomialRing(Qp);
    rQp:=Qpx!r;
    zero:=HenselLift(rQp,x1);
    sQp:=rQp div (Qpx.1-zero);
    denom:=Qt!Kt!((Qt!OKt!(xt-Coefficient(xt,0)))^(-1)*(Qt!Kt!(1/Evaluate(sQp,xt))));
  end if;

  // determine the number of points at infinity

  FF:=fun_field(data);
  infplaces:=InfinitePlaces(FF);
  infpoints:=0;
  for i:=1 to #infplaces do
    infpoints:=infpoints+Degree(infplaces[i]);
  end for;

  // compute the powerseries expansions of the basis elements of H^1(Y)

  derxt:=Qt!Kt!Derivative(xt); 
  diffs:=[];
  for i:=1 to 2*g+infpoints-1 do
    basisxt:=Evaluate(basis[i],xt);
    for j:=1 to d do
      basisxt[1][j]:=Qt!Kt!basisxt[1][j];
    end for;
    diffs[i]:=InnerProduct(Vector(basisxt*derxt*lc_r*denom),bt);
    diffs[i]:=Qt!Kt!diffs[i];
    if Coefficient(diffs[i],-1) ne 0 then
      diffs[i]:=diffs[i]-Coefficient(diffs[i],-1)*Kt.1^(-1);
    end if;
  end for;

  omegax:=[];
  for i:=1 to 2*g+infpoints-1 do
    omegax[i]:=diffs[i]/denombasis; // scaling by denominator of our basis
  end for;

  // tiny single integrals

  Omegax:=[];
  for i:=1 to 2*g do
    Omegax[i]:=-Integral(omegax[i]);
  end for;

  singleintegrals:=[];
  for i:=1 to 2*g do
    if index eq 0 then // x-x(P1) is the local coordinate
      singleintegrals[i]:=Evaluate(Omegax[i],x2-x1);
    else // b[index]-b[index](P1) is the local coordinate
      singleintegrals[i]:=Evaluate(Omegax[i],b2[index]-b1[index]);
    end if;
  end for;
 
  // tiny double integral

  dgx:=0;
  for i:=1 to 2*g do
    for j:=1 to 2*g do
      dgx:=dgx+omegax[i]*(K!Z[i,j])*Omegax[j]; 
    end for;
  end for;
  for i:=1 to infpoints-1 do
    dgx:=dgx-eta[i]*omegax[2*g+i];
  end for;  
  gx:=Integral(dgx);

  if index eq 0 then // x-x(P1) is the local coordinate
    doubleintegral:=Evaluate(gx,x2-x1);
  else // b[index]-b[index](P1) is the local coordinate
    doubleintegral:=Evaluate(gx,b2[index]-b1[index]);
  end if;

  C:=IdentityMatrix(K,2*g+2);
  // entries in first column (except last one) 
  for i:=1 to 2*g do
    C[i+1,1]:=-singleintegrals[i];
  end for;
  // entries in the last row (except first one)
  for i:=1 to 2*g do
    for j:=1 to 2*g do
      C[2*g+2,i+1]:=C[2*g+2,i+1]-singleintegrals[j]*(K!Z[j,i]);
    end for;
  end for;
  // bottom left entry
  C[2*g+2,1]:=-doubleintegral;

  return C;

end function;



