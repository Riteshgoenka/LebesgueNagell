D := 3;
_<x> := PolynomialRing(Rationals());
F<r> := NumberField(x^2 - D);
OF := Integers(F);
M := HilbertCuspForms(F, ((r-1)*OF)^10*(r*OF)^2);
Mnew := NewSubspace(M);
print(Dimension(MNew));
N := NewformDecomposition(Mnew);
print(#N);

// This function computes the bound \mathcal{B}_q(f).
function Bound(f,q)
  b := 1;
  B := 1;
  factQ:=Factorisation(q*OF);
  for a in [1..q] do
    if [a,b] ne [q,q] then
      L:=[];
      C:= EllipticCurve([0,4*(4*r)*b,0,2*(4*r)*(a+(4*r)*b^2),0]);
      for i in [1..#factQ] do
        Q:=factQ[i,1];
        if LocalInformation(C,Q)[3] eq 0 then
          L:=Append(L,Integers()!Norm(TraceOfFrobenius(C,Q)-HeckeEigenvalue(f,Q)));
        else
          L:=Append(L,Integers()!Norm((Norm(Q)+1)^2-HeckeEigenvalue(f,Q)^2));
        end if;
      end for;
      B:=B*Gcd(L);
    end if;
  end for;
  return q*B;
end function;

S := {};
for i in [1..#N] do
  V := [];
  g := Eigenform(N[i]);
  for q in PrimesInInterval(5,18) do
    Include(~V,Bound(g,q));
  end for;
  h := Gcd(V);
  if h eq 0 then
    print i, "Found it!";
  else
    Fac := Factorisation(h);
    for pr in Fac do
      Include(~S,pr[1]);
    end for;
    print Fac;
  end if;
end for;
print S;