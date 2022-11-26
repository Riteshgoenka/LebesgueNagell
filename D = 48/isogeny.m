K<z>:=QuadraticField(3);
S<X>:=PolynomialRing(K);

L<z1>:=ext<K|X^2-2>;
r:=L!(4 * z);

R:=PolynomialRing(L,2);
F<a,b>:=FieldOfFractions(R);

/* Frey curve for a^2 - 48 b^4 = c^p. */
E:=EllipticCurve([0,4*r*b,0,2*r*(a+r*b^2),0]);

/* Conjugate curve */
Ec:=EllipticCurve([0,-4*r*b,0,-2*r*(a-r*b^2),0]);

/* There is a cyclic subgroup of order 2 defined over L (in fact over K) on Ec. */
f2:=DivisionPolynomial(Ec,2);
f2l:=Factorization(f2);
E2,phi2:=IsogenyFromKernel(Ec,f2l[1][1]);

/* Unfortunately, E2 and E are not isomorphic over K = Q(3^(1/2)), but only over L = Q(2^(1/2),3^(1/2)). */
IsIsomorphic(E2,E);
