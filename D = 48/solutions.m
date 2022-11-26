D := 3;
_<x> := PolynomialRing(Rationals());
F<r> := NumberField(x^2 - D);
OF := Integers(F);
M := HilbertCuspForms(F, ((r-1)*OF)^10*(r*OF)^2);
Mnew := NewSubspace(M);
N := NewformDecomposition(Mnew);
F1 := Eigenform(N[41]);
F2 := Eigenform(N[56]);

a := -7;
b := 1;
C := EllipticCurve([0,4*(4*r)*b,0,2*(4*r)*(a+(4*r)*b^2),0]);
for q in PrimesInInterval(5,30) do
	factQ := Factorisation(q*OF);
	for i in [1..#factQ] do
        Q := factQ[i,1];
		if LocalInformation(C,Q)[3] eq 0 then
			print TraceOfFrobenius(C,Q)-HeckeEigenvalue(F2,Q);
		else
    		print (Norm(Q)+1)^2-HeckeEigenvalue(F2,Q)^2;
		end if;
	end for;
end for;

a := 7;
b := 1;
C := EllipticCurve([0,4*(4*r)*b,0,2*(4*r)*(a+(4*r)*b^2),0]);
for q in PrimesInInterval(5,30) do
	factQ := Factorisation(q*OF);
	for i in [1..#factQ] do
        Q := factQ[i,1];
		if LocalInformation(C,Q)[3] eq 0 then
			print TraceOfFrobenius(C,Q)-HeckeEigenvalue(F1,Q);
		else
    		print (Norm(Q)+1)^2-HeckeEigenvalue(F1,Q)^2;
		end if;
	end for;
end for;

a := -7;
b := 1;
C := EllipticCurve([0,4*(4*r)*b,0,2*(4*r)*(a+(4*r)*b^2),0]);
for q in PrimesInInterval(5,30) do
	factQ := Factorisation(q*OF);
	for i in [1..#factQ] do
        Q := factQ[i,1];
		if LocalInformation(C,Q)[3] eq 0 then
			print TraceOfFrobenius(C,Q)-HeckeEigenvalue(F1,Q);
		else
    		print (Norm(Q)+1)^2-HeckeEigenvalue(F1,Q)^2;
		end if;
	end for;
end for;

a := 7;
b := 1;
C := EllipticCurve([0,4*(4*r)*b,0,2*(4*r)*(a+(4*r)*b^2),0]);
for q in PrimesInInterval(5,30) do
	factQ := Factorisation(q*OF);
	for i in [1..#factQ] do
        Q := factQ[i,1];
		if LocalInformation(C,Q)[3] eq 0 then
			print TraceOfFrobenius(C,Q)-HeckeEigenvalue(F2,Q);
		else
    		print (Norm(Q)+1)^2-HeckeEigenvalue(F2,Q)^2;
		end if;
	end for;
end for;
