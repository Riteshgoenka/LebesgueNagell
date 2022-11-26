
/* Computing conductor at 2. */

/* We use Magma's routine to compute conductor. */

R<X>:=PolynomialRing(RationalField());

K<z>:=NumberField(X^2-80);

/* Prime above (2). */
O:=MaximalOrder(K);
L:=Factorization(2*O);
P:=L[1][1];
pi:=UniformizingElement(P);
k,f:=ResidueClassField(P);

p:=Place(P);
R<v,t>:=FieldOfFractions(PolynomialRing(K,2));
S<x,y>:=PolynomialRing(R,2);

M:=2^8;

print "Computing exponent of conductor at prime above 2 for E.";
print "Going through residue classes of (a,b) modulo",M;
print "Let q_2 be the prime above 2 in K_beta.";
print "Let f be the exponent of the conductor at q_2.";
print "Let c be the Kodaira symbol at q_2.";
print "Let d be the valuation of the discriminant at q_2 of initial model.";
print "Let dm be the valuation of the discriminant at q_2 of the minimal model.";

printf "%3o %3o %3o %3o %3o %3o\n", "a","b","f", "c","d","dm";

RM:=ResidueClassRing(M);
listk:=[[1,1]];

for a:=0 to M-1 do
  for b:=1 to 1 do
    if ([a mod 2, b mod 2] in listk) then
      E:=EllipticCurve([0,0,0,-3*z*(4*a+5*z),2*z*(2*a^2+11*z^2+14*z*a)]);
      LI:=LocalInformation(E,P);
      printf "%3o %3o %3o %3o %3o %3o\n",a,b,LI[3],LI[5],Valuation(Discriminant(E),p),LI[2];
    end if;
  end for;
end for;


