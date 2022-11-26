
/* Computing conductor at 3. */

/* We use Magma's routine to compute conductor. */

R<X>:=PolynomialRing(RationalField());

K<z>:=NumberField(X^2-48);

/* Prime above (3). */
O:=MaximalOrder(K);
L:=Factorization(3*O);
P:=L[1][1];
pi:=UniformizingElement(P);
k,f:=ResidueClassField(P);

p:=Place(P);
R<v,t>:=FieldOfFractions(PolynomialRing(K,2));
S<x,y>:=PolynomialRing(R,2);

M:=3^4;

print "Computing exponent of conductor at prime above 3 for E.";
print "Going through residue classes of (a,b) modulo",M;
print "Let q_3 be the prime above 3 in K_beta.";
print "Let f be the exponent of the conductor at q_3.";
print "Let c be the Kodaira symbol at q_3.";
print "Let d be the valuation of the discriminant at q_3 of initial model.";
print "Let dm be the valuation of the discriminant at q_3 of the minimal model.";

printf "%3o %3o %3o %3o %3o %3o\n", "a","b","f", "c","d","dm";

RM:=ResidueClassRing(M);
listk:=[[0,1],[1,1],[2,1]];

for a:=0 to M-1 do
  for b:=1 to 1 do
    if ([a mod 3, b mod 3] in listk) then
      E:=EllipticCurve([0,4*z,0,2*z*(a+z),0]);
      LI:=LocalInformation(E,P);
      printf "%3o %3o %3o %3o %3o %3o\n",a,b,LI[3],LI[5],Valuation(Discriminant(E),p),LI[2];
    end if;
  end for;
end for;


