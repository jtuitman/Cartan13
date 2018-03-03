//////////////////////////////////////////////
// Functions for computing Hodge structures //
//////////////////////////////////////////////


find_eta:=function(data,denombasis,Z)

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

  infplacesKinf:=InfinitePlaces(FFKinf); // places at infinity all of degree 1

  L:=[];
  for i:=1 to (2*g+#infplacesKinf-1) do
    fun:=FFKinf!0;
    for j:=1 to d do
      fun:=fun+Evaluate(basis[i][j],Kinfx.1)*b0fun[j];
    end for;
    fun:=fun/denombasis;
    L[i]:=fun;
  end for;

  // set up the linear system eta*A=v

  v:=[];
  A:=ZeroMatrix(Kinf,#infplacesKinf-1,#infplacesKinf);

  for i:=1 to #infplacesKinf do
    
    P:=infplacesKinf[i];
    dxdt:=Derivative(Expand(FFKinf!Kinfx.1,P));
    zinv:=Expand(LeadingCoefficient(r)/(FFKinf!Evaluate(r,Kinfx.1)),P);
    omega:=[];
    Omega:=[];
    for j:=1 to 2*g do
      omega[j]:=Expand(L[j],P)*dxdt*zinv; // expansion of omega_j at P
      Omega[j]:=Integral(omega[j]);       // expansion of Omega_j at P
    end for;
    
    omegaZOmega:=0; // omega*Z*Omega
    for i:=1 to 2*g do
      for j:=1 to 2*g do
        omegaZOmega:=omegaZOmega+omega[i]*Z[i,j]*Omega[j];
      end for;
    end for;

    OmegaZomega:=0; // Omega*Z*omega
    for i:=1 to 2*g do
      for j:=1 to 2*g do
        OmegaZomega:=OmegaZomega+Omega[i]*Z[i,j]*omega[j];
      end for;
    end for;
 
    OmegaZTomega:=0; // Omega*Tranpose(Z)*omega
    for i:=1 to 2*g do
      for j:=1 to 2*g do
        OmegaZTomega:=OmegaZTomega+Omega[i]*Z[j,i]*omega[j];
      end for;
    end for;

    v[i]:=Coefficient(omegaZOmega-OmegaZomega+OmegaZTomega,-1); // residue of eta at P

    for j:=1 to (#infplacesKinf-1) do
      omega[2*g+j]:=Expand(L[2*g+j],P)*dxdt*zinv;
      A[j,i]:=Coefficient(omega[2*g+j],-1); // residue of omega_{2g+j} at P
    end for;

  end for;

  eta:=Solution(A,Vector(v)); // solve for eta

  return eta;

end function;
