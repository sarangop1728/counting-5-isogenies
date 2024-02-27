
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

poly := 5*psi5facts[1][1]; // quadratic factor of psi5 in ZZ[T]

// 5*T^2 - 300*c*T - 18954*a - 15750*c^2 + 7128*b

cs := Coefficients(IsogenyFromKernel(E, poly)); // 5-isogenous curve

// The GCD of the coefficients of cs[4] is 33750 = 2*3*(3*25)^2.
// The GCD of the coefficients of cs[6] is 3375000 = 2^3*(3*25)^3;

u := 1/(3*25);

EE := EllipticCurve([cs[4]*u^2, cs[5]*u^3]); // simplified equation 

// Elliptic Curve defined by 
// y^2 = x^3 + (-18*a - 30*c^2 + 36*b)*x + (-144*c*a - 160*c^3 + 72*b*c)
// over Ring of Fractions of Affine Algebra of rank 3 over Rational Field