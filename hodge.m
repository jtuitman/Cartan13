//////////////////////////////////////////////
// Functions for computing Hodge structures //
//////////////////////////////////////////////

hodge_data:=function(data,denombasis,Z,bpt)

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

  b0funKinf:=[]; // functions b^0 (finite integral basis)
  for i:=1 to d do
    b0i:=FFKinf!0;
    for j:=1 to d do
      b0i:=b0i+Evaluate(W0[i,j],Kinfx.1)*FFKinf.1^(j-1);
    end for;
    b0funKinf[i]:=b0i;
  end for;

  infplacesKinf:=InfinitePlaces(FFKinf); // places at infinity all of degree 1, will be denoted by P

  L:=[];
  for i:=1 to (2*g+#infplacesKinf-1) do
    fun:=FFKinf!0;
    for j:=1 to d do
      fun:=fun+Evaluate(basis[i][j],Kinfx.1)*b0funKinf[j];
    end for;
    fun:=fun/denombasis;
    L[i]:=fun;
  end for;

  // compute the expansions omega_x, Omega_x, b^0_x

  omegax:=[]; // expansions of omega
  Omegax:=[]; // expansions of Omega
  b0funx:=[]; // expansions of b^0
  xfunx:=[];  // expansions of x

  for i:=1 to #infplacesKinf do

    P:=infplacesKinf[i];
    
    xfunx[i]:=Expand(FFKinf!Kinfx.1,P);
    dxdt:=Derivative(xfunx[i]);
    
    zinv:=Expand(LeadingCoefficient(r)/(FFKinf!Evaluate(r,Kinfx.1)),P);
    
    omegaP:=[];
    for j:=1 to 2*g+#infplacesKinf-1 do
      omegaP[j]:=Expand(L[j],P)*dxdt*zinv;
    end for;
    omegax:=Append(omegax,omegaP);
    
    OmegaP:=[];
    for j:=1 to 2*g do
      OmegaP[j]:=Integral(omegaP[j]);
    end for;
    Omegax:=Append(Omegax,OmegaP);

    b0funP:=[];
    for j:=1 to d do
      b0funP[j]:=Expand(b0funKinf[j],P);
    end for;
    b0funx:=Append(b0funx,b0funP);

  end for;

  // compute expansions of Omega*Z*omega at all points at infinity

  OmegaZomega:=[];
  for i:=1 to #infplacesKinf do
    OmegaZomegaP:=0;
    for j:=1 to 2*g do
      for k:=1 to 2*g do
        OmegaZomegaP:=OmegaZomegaP+omegax[i][j]*Z[j,k]*Omegax[i][k];
      end for;
    end for;
    OmegaZomega:=Append(OmegaZomega,OmegaZomegaP);
  end for;

  // set up the linear system eta*A=v satisfied by eta
  
  v:=[];
  A:=ZeroMatrix(Kinf,#infplacesKinf-1,#infplacesKinf);
  for i:=1 to #infplacesKinf do
    v[i]:=-Coefficient(OmegaZomega[i],-1); // residue of eta at i-th point of infinity
    for j:=1 to #infplacesKinf-1 do
      A[j,i]:=Coefficient(omegax[i][2*g+j],-1); // residue of omega_{2g+j} at i-th point at infinity
    end for;
  end for;

  eta:=Solution(A,Vector(v)); // solve for eta

  gx:=[]; // functions g_x
  for i:=1 to #infplacesKinf do
    dgxi:=OmegaZomega[i]; 
    for j:=1 to (#infplacesKinf-1) do
      dgxi:=dgxi+eta[j]*omegax[i][2*g+j]; 
    end for;
    gx[i]:=Integral(dgxi);
  end for;

  poleorder:=0;
  for i:=1 to #infplacesKinf do
    val:=Valuation(gx[i]);
    for j:=1 to 2*g do
      val:=Minimum(val,Valuation(Omegax[i][j]));
    end for;
    poleorder:=Minimum(poleorder,val);
  end for;

  done:=false;
  degx:=0;
  while not done do // try larger and larger degree in x
    
    for i:=1 to #infplacesKinf do
      for j:=1 to d do
        poleorder:=Minimum(poleorder,Valuation(b0funx[i][j])+degx*Valuation(xfunx[i]));
      end for;
    end for;

    v:=[]; // coefficients of principal parts of all gx 
    cnt:=0;
    for i:=1 to #infplacesKinf do
      for j:=poleorder to -1 do
        cnt:=cnt+1;
        v[cnt]:=Coefficient(gx[i],j);
      end for;
    end for;

    rows:=[];

    for i:=1 to g do
      row:=[];
      cnt:=0;
      for j:=1 to #infplacesKinf do
        for k:=poleorder to -1 do
          cnt:=cnt+1;
          row[cnt]:=Coefficient(Omegax[j][i+g],k); // coefficients of principal part of Omegax_{i+g} at jth point at infinity
        end for; 
      end for;
      rows:=Append(rows,row);
    end for;

    for i:=1 to d do
      for j:=0 to degx do
        row:=[];
        cnt:=0;
        for k:=1 to #infplacesKinf do
          for l:=poleorder to -1 do
            cnt:=cnt+1;
            row[cnt]:=Coefficient(b0funx[k][i]*xfunx[k]^j,l); // coefficients of principal part of x^j*b^0_i at kth point at infinity
          end for;
        end for;
        rows:=Append(rows,row);  
      end for;
    end for;   
      
    suc,sol:=IsConsistent(Matrix(rows),Vector(v));
    if suc then
      done:=true;
    else // if no success, increase the degree in x
      degx:=degx+1;
    end if;
  
  end while;

  // read off beta from solution

  beta:=[];
  for i:=1 to g do
    beta[i]:=sol[i];
  end for;

  // read off gamma from solution

  Qt:=PolynomialRing(RationalField());
  gamma:=[];
  cnt:=g;
  for i:=1 to d do
    poly:=Qx!0;
    for j:=0 to degx do
      cnt:=cnt+1;
      poly:=poly+(RationalField()!sol[cnt])*Qt.1^j;
    end for;
    gamma:=Append(gamma,poly);
  end for;

  b0fun:=[]; // functions b^0 (finite integral basis)
  for i:=1 to d do
    b0i:=FF!0;
    for j:=1 to d do
      b0i:=b0i+Evaluate(W0[i,j],Qx.1)*FF.1^(j-1);
    end for;
    b0fun[i]:=b0i;
  end for;

  // substract constant such that gamma(bpt)=0

  gamma_FF:=FF!0;
  for i:=1 to d do
    gamma_FF:=gamma_FF+Evaluate(gamma[i],Qx.1)*b0fun[i];
  end for;
  gamma[1]:=gamma[1]-Evaluate(gamma_FF,bpt); 

  // TODO analyse t-adic precision
  // TODO beta,gamma off by factor 3/2 compared to paper

  return eta,beta,gamma;

end function;
