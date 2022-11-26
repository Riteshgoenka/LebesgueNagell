K<z2p>:=QuadraticField(5);
S<X>:=PolynomialRing(K);

L<z>:=ext<K|6*X^2-5-z2p>;
r:=L!(4 * z2p);

gamma:=-18*z^3 + 15*z^2 + 10*z - 5;
gammac:=-24*z^3 - 15*z^2 + 30*z + 20;

R:=PolynomialRing(L,2);
F<a,b>:=FieldOfFractions(R);

/* Frey curve for a^2 - 2 b^6 = c^p. */
E:=EllipticCurve([0,0,0,-3*r*(4*a+5*r*b^3)*b*gamma^2,2*r*(2*a^2+11*r^2*b^6+14*r*a*b^3)*gamma^3]);

/* Conjugate curve */
Ec:=EllipticCurve([0,0,0,3*r*(4*a-5*r*b^3)*b*gammac^2,-2*r*(2*a^2+11*r^2*b^6-14*r*a*b^3)*gammac^3]);

/* There is a cyclic subgroup of order 3 over L on Ec. */
f3:=DivisionPolynomial(Ec,3);
f3l:=Factorization(f3);
P3:=f3l[1][1];

/* Find using Velu the elliptic curve E3 = E/P3. */
E3,phi3:=IsogenyFromKernel(Ec,P3);

/* E3 and E are isomorphic over L = Q(((5+5^(1/2))/6)^(1/2)). */
IsIsomorphic(E3,E);

