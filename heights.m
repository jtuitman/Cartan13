/////////////////////////////////////
// functions for computing heights //
/////////////////////////////////////


function frob_equiv_iso(G,data)

  // given a matrix G defining Frob^-1 on a mixed extension this returns s^phi

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
