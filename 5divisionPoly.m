
// This code gives the quadratic factor of
// the 5-division polynomial corresponding
// to (a,b,c) satisfying a^2 + b^2 = c^4.

S<x,y,z> := PolynomialRing(Rationals(), 3);
R := quo<S | x^2 + y^2 - z^4>;
FF<a,b,c> := FieldOfFractions(R);

e4 := 246*a + 228*b + 250*c^2;
e6 := -c*(10008*a + 1044*b + 10000*c^2);

E := EllipticCurve([-27*e4,-54*e6]);
psi5 := DivisionPolynomial(E,5);
psi5facts := Factorization(psi5);

psi5facts[1];
// T^2 - 60*T - 18954/5*a + (-15750*c^2 + 7128*b)/5