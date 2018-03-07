//////////////////////////////////////////////
// Functions for computing Hodge structures //
//////////////////////////////////////////////


hodge_data:=function(data,denombasis,Z)

  // Compute the 1-form eta, as a vector of coefficients
  // w.r.t. basis[i]/denombasis for i=2g+1,...,2g+k-1 where
  // k is the number of points lying over x=infinity.

  Q:=data`Q; g:=data`g; r:=data`r; W0:=data`W0; basis:=data`basis;

  d:=Degree(Q);

  // find the points at infinity:

  Qx:=RationalFunctionField(RationalField()); Qxy:=PolynomialRing(Qx);

  f:=Qxy!0;
  for i:=0 to d do
    for j:=0 to Degree(Coefficient(Q,i)) do
      f:=f+Coefficient(Coefficient(Q,i),j)*Qxy.1^i*Qx.1^j;
    end for;
  end for;  
  FF:=FunctionField(f); // function field of curve

  infplaces:=InfinitePlaces(FF);

  Kinf:=RationalField();
  for i:=1 to #infplaces do
    Kinf:=Compositum(Kinf,ResidueClassField(infplaces[i]));
  end for;

  Kinfx:=RationalFunctionField(Kinf); Kinfxy:=PolynomialRing(Kinfx);

  finf:=Kinfxy!0;
  for i:=0 to d do
    for j:=0 to Degree(Coefficient(Q,i)) do
      finf:=finf+Coefficient(Coefficient(Q,i),j)*Kinfxy.1^i*Kinfx.1^j;
    end for;
  end for;  
  FFKinf:=FunctionField(finf); // function field of curve over Kinf

  b0fun:=[]; // functions b^0 (finite integral basis)
  for i:=1 to d do
    b0i:=FFKinf!0;
    for j:=1 to d do
      b0i:=b0i+Evaluate(W0[i,j],Kinfx.1)*FFKinf.1^(j-1);
    end for;
    b0fun[i]:=b0i;
  end for;

  infplacesKinf:=InfinitePlaces(FFKinf); // places at infinity all of degree 1, will be denoted by P

  L:=[];
  for i:=1 to (2*g+#infplacesKinf-1) do
    fun:=FFKinf!0;
    for j:=1 to d do
      fun:=fun+Evaluate(basis[i][j],Kinfx.1)*b0fun[j];
    end for;
    fun:=fun/denombasis;
    L[i]:=fun;
  end for;

  // compute the expansions omega_x and Omega_x

  omegax:=[];
  Omegax:=[];
  for i:=1 to #infplacesKinf do
    P:=infplacesKinf[i];
    dxdt:=Derivative(Expand(FFKinf!Kinfx.1,P));
    zinv:=Expand(LeadingCoefficient(r)/(FFKinf!Evaluate(r,Kinfx.1)),P);
    omegaP:=[];
    OmegaP:=[];
    for j:=1 to 2*g+#infplacesKinf-1 do
      omegaP[j]:=Expand(L[j],P)*dxdt*zinv;
    end for;
    omegax:=Append(omegax,omegaP);
    for j:=1 to 2*g do
      OmegaP[j]:=Integral(omegaP[j]);
    end for;
    Omegax:=Append(Omegax,OmegaP);
  end for;

  // set up the linear system eta*A=v
  
  v:=[];
  A:=ZeroMatrix(Kinf,#infplacesKinf-1,#infplacesKinf);

  omegaZOmega_minus_OmegaZminusZTomega:=[]; // expansions of omega*Z*Omega-Omega*(Z-Z^T)*omega at all P

  for i:=1 to #infplacesKinf do
    
    omegaZOmega:=0; // omega*Z*Omega
    for j:=1 to 2*g do
      for k:=1 to 2*g do
        omegaZOmega:=omegaZOmega+omegax[i][j]*Z[j,k]*Omegax[i][k];
      end for;
    end for;

    OmegaZminusZTomega:=0; // Omega*(Z-Z^T)*omega
    for j:=1 to 2*g do
      for k:=1 to 2*g do
        OmegaZminusZTomega:=OmegaZminusZTomega+Omegax[i][j]*(Z[j,k]-Z[k,j])*omegax[i][k];
      end for;
    end for;

    omegaZOmega_minus_OmegaZminusZTomega[i]:=omegaZOmega-OmegaZminusZTomega; // expansion of omega*Z*Omega-Omega*(Z-Z^T)*omega at i-th point at infinity

    v[i]:=Coefficient(omegaZOmega_minus_OmegaZminusZTomega[i],-1); // TODO clear up sign
    for j:=1 to #infplacesKinf-1 do
      A[j,i]:=Coefficient(omegax[i][2*g+j],-1); // residue of omega_{2g+j} at P
    end for;

  end for;

  eta:=-Solution(A,Vector(v)); // solve for eta, TODO clear up sign

  gx:=[]; // functions g_x
  for i:=1 to #infplacesKinf do
    dgxi:=omegaZOmega_minus_OmegaZminusZTomega[i]; 
    for j:=1 to (#infplacesKinf-1) do
      dgxi:=dgxi+eta[j]*omegax[i][2*g+j]; // TODO clear up sign
    end for;
    gx[i]:=Integral(dgxi);
  end for;

  // TODO gammaFil and betaFil

  return eta,gx;

end function;
