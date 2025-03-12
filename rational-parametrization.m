// We compute a rational parametrization
// for the modular curve X0(5).

_<t> := FunctionField(Rationals());
// Taken from pg. 1247 of Halberstadt:
j := (t^2 + 10*t + 5)^3/t;
E := EllipticCurveWithjInvariant(j);
E := SimplifiedModel(E);
cs := Coefficients(E);
Factorization(Numerator(cs[4]));
// [
//     <$.1^2 + 10*$.1 + 5, 3>
// ]
Factorization(Denominator(cs[4]));
// [
//     <$.1^2 + 4*$.1 - 1, 2>,
//     <$.1^2 + 22*$.1 + 125, 1>
// ]
Factorization(Numerator(cs[5]));
// [
//     <$.1^2 + 10*$.1 + 5, 3>
// ]
Factorization(Denominator(cs[5]));
// [
//     <$.1^2 + 4*$.1 - 1, 2>,
//     <$.1^2 + 22*$.1 + 125, 1>
// ]

u := (t^2 + 4*t -1)*(t^2 + 22*t + 125)/(t^2 + 10*t + 5);
E := EllipticCurve([cs[4]*u^2, cs[5]*u^3]);
P<T> := PolynomialRing(Rationals());
cs := Coefficients(E);
f0 := P!cs[4];
g0 := P!cs[5];
lcm_den := LeastCommonMultiple([Denominator(c) : c in Coefficients(f0) cat Coefficient\
s(g0)]);
Factorization(lcm_den);
// [ <2, 5>, <3, 3> ]
cs := Coefficients(E);
v := 2^2*3;
E := EllipticCurve([cs[4]*v^2, cs[5]*v^3]); // has coefficients in Q[t]
cs := Coefficients(E);
f := P!cs[4];
Factorization(f);
// f = -3(T^2 + 10*T + 5)(T^2 + 22*T + 125) : degree 4

g := P!cs[5];
Factorization(g);
// g = 2(T^2 + 4*T - 1)(T^2 + 22*T + 125)^2 : degree 6

// Consistency check:
// The 5-division polynomial factors.
phi5facts := Factorization(DivisionPolynomial(E,5));
phi5facts[1];
