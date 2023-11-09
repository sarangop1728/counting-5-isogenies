load "canring.m";

// EQUATIONS FOR e4 AND e6 in terms of a,b,c satisfying a^2 + b^2 = c^4:

R, I, MGamma, RmodI, genlist, genweights := CanonicalRing(5);
MonomialsOfWeightedDegree(R, 4);
_<q> := Parent(genlist[1]);
f2 := qExpansion(genlist[1],100); // degree 2
f4 := qExpansion(genlist[2],100); // degree 4
g4 := qExpansion(genlist[3],100); // degree 4
e4 := Eisenstein(4,q : Precision := 100);


a := 2*g4 - 1/4*f2^2; 
b := f4 + 7*g4; 
c := 1/2*f2;
// check that a^2 + b^2 = c^4:
MonomialsOfWeightedDegree(R, 8);
Kernel(Transpose(Matrix([[Coefficient(g,n) : g in [a^2,b^2,c^4]]: n in [0..20]])));
a^2 + b^2 eq c^4;

// EQUATIONS FOR e4 and e6 with respect to a,b,c
MonomialsOfWeightedDegree(R,4);
Kernel(Transpose(Matrix([[Coefficient(g,n) : g in [e4,a,b,c^2]]: n in [0..20]])));
e4 eq 246*a + 228*b + 250*c^2;
e6 := Eisenstein(6,q : Precision := 100);
MonomialsOfWeightedDegree(R,6);
Kernel(Transpose(Matrix([[Coefficient(g,n) : g in [e6,a*c,b*c,c^3]]: n in [0..20]])));
e6 eq -10008*a*c - 1044*b*c - 10000*c^3;

// Universal elliptic curve over M0(5):

_<t> := FunctionField(Rationals());
j := (t^2 + 10*t + 5)^3/t; // taken from pg. 1247 of Halberstadt
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

u := (t^2 + 4*t -1)*(t^2 + 22*t + 125);
E := EllipticCurve([cs[4]*u^2, cs[5]*u^3]); // has coefficients in Q[t]
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
// f =
// [
//     <T^2 + 10*T + 5, 3>,
//     <T^2 + 22*T + 125, 1>
// ]
g := P!cs[5];
Factorization(g);
// g =
// [
//     <T^2 + 4*T - 1, 1>,
//     <T^2 + 10*T + 5, 3>,
//     <T^2 + 22*T + 125, 2>
// ]
phi5facts := Factorization(DivisionPolynomial(E,5));
phi5facts[1];