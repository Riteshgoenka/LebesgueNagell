N:=3^3*5^2;
 
G3<a3>:=FullDirichletGroup(3);
G5<a5>:=FullDirichletGroup(5);

G:=DirichletGroup(N, CyclotomicField(4));

eps:=G!(a3*a5);

SetVerbose("ModularSymbols",0);
FullSym := ModularSymbols(eps^(-1),2);
CuspSym := CuspidalSubspace(FullSym);
NewSym := NewSubspace(CuspSym);
NewList := NewformDecomposition(NewSym);
SortList := SortDecomposition(NewList);

datafilename := "all-675.data";

B:=5000;

SetVerbose("ModularSymbols",0);

F:=Open(datafilename,"w");

fprintf F,"// qexplist is a list of items [[m],[qexp]] \n";
fprintf F,"// m is the modulus for the field of coefficients of the form \n";
fprintf F,"// qexp is the sequence of a_p over p prime for the form \n";
fprintf F,"// p runs through primes <= %o \n \n",B;

fprintf F,"R<x>:=PolynomialRing(CyclotomicField(4)); \n";
fprintf F,"qexplist:=[ \n";

lk:=[];

for k:=1 to #SortList do
  f:=AssociatedNewSpace(SortList[k]);
  qexp:=SystemOfEigenvalues(f,B);

  S:=Parent(qexp[1]);
  R<x>:=PolynomialRing(CyclotomicField(4));
  if S cmpeq RationalField() then
    m:="x";
    fprintf F, "[%o,%o], \n",[m],qexp;
  else
    m:=DefiningPolynomial(S);
	print(m);
    if (k ne #SortList) then
      fprintf F, "[%o,%o], \n",[m],[R!ElementToSequence(aq):aq in qexp];
    else
      fprintf F, "[%o,%o] \n",[m],[R!ElementToSequence(aq):aq in qexp];
    end if;
  end if;
end for;

fprintf F,"]; \n";

delete(F);