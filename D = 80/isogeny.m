K<z>:=QuadraticField(5);
S<X>:=PolynomialRing(K);

L<z1>:=ext<K|X^2-3>;
r:=L!(4 * z);

R:=PolynomialRing(L,2);
F<a,b>:=FieldOfFractions(R);

/* Frey curve for a^2 - 2 b^6 = c^p. */
E:=EllipticCurve([0,0,0,-3*r*(4*a+5*r*b^3)*b,2*r*(2*a^2+11*r^2*b^6+14*r*a*b^3)]);

/* Conjugate curve */
Ec:=EllipticCurve([0,0,0,3*r*(4*a-5*r*b^3)*b,-2*r*(2*a^2+11*r^2*b^6-14*r*a*b^3)]);

/* There is a cyclic subgroup of order 3 defined over L (in fact over K) on Ec. */
f3:=DivisionPolynomial(Ec,3);
f3l:=Factorization(f3);
E3,phi3:=IsogenyFromKernel(Ec,f3l[1][1]);

/* Unfortunately, E3 and E are not isomorphic over K = Q(5^(1/2)), but only over L = Q(3^(1/2),5^(1/2)). */
IsIsomorphic(E3,E);
