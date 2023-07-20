
load "canring.m";

// EQUATIONS FOR e4 AND e6:

R, I, MGamma, RmodI, genlist, genweights := CanonicalRing(5);
R<cp,bp,ap> := R;  // names
MonomialsOfWeightedDegree(R, 4);
_<q> := Parent(genlist[1]);
f2 := qExpansion(genlist[1],100); // cp
f4 := qExpansion(genlist[2],100); // bp?
g4 := qExpansion(genlist[3],100); // ap?
e4 := Eisenstein(4,q : Precision := 100);
Kernel(Transpose(Matrix([[Coefficient(g,n) : g in [f2^2,f4,g4,e4]] : n in [0..20]])));
e4 eq f2^2+228*f4+2088*g4; // equation for e4
e6 := Eisenstein(6,q : Precision := 100);
MonomialsOfWeightedDegree(R, 6);
Kernel(Transpose(Matrix([[Coefficient(g,n) : g in [f2^3,f2*f4,f2*g4,e6]] : n in [0..20]])));
e6 eq f2^3-522*f2*f4-13662*f2*g4; // equation for e6


// UNIVERSAL ELLIPTIC CURVE OVER M0(5):

_<t> := FunctionField(Rationals());
j := 1/2*(t^6*(19*t^2+8*t-16)^3/(t^11*(3*t-2)));  // from LMFDB
E := EllipticCurveWithjInvariant(j);
E := SimplifiedModel(E);
cs := Coefficients(E);
Factorization(Numerator(cs[4]));
// [
//     <$.1^2 + 8/19*$.1 - 16/19, 3>
// ]
Factorization(Denominator(cs[4]));
// [
//     <$.1^2 - 32/11*$.1 + 16/11, 2>,
//     <$.1^2 + 40/29*$.1 + 16/29, 1>
// ]
Factorization(Numerator(cs[5]));
// [
//     <$.1^2 + 8/19*$.1 - 16/19, 3>
// ]
Factorization(Denominator(cs[5]));
// [
//     <$.1^2 - 32/11*$.1 + 16/11, 2>,
//     <$.1^2 + 40/29*$.1 + 16/29, 1>
// ]
u := (t^2-32/11*t+16/11)*(t^2+40/29*t+16/29)/(t^2+8/19*t-16/19);
E := EllipticCurve([cs[4]*u^2, cs[5]*u^3]);
cs := Coefficients(E);
Factorization(Numerator(cs[4]));
// [
//     <$.1^2 + 8/19*$.1 - 16/19, 1>,
//     <$.1^2 + 40/29*$.1 + 16/29, 1>
// ]
Factorization(Numerator(cs[5]));
// [
//     <$.1^2 - 32/11*$.1 + 16/11, 1>,
//     <$.1^2 + 40/29*$.1 + 16/29, 2>
// ]
_<a,b,c> := PolynomialRing(Rationals(),3);
Factorization(Parent(b)!(Evaluate(Coefficients(E)[4],b/(4*(c-2)))*(c-2)^4)); 
// [
//     <b^2 + 32/19*b*c - 64/19*b - 256/19*c^2 + 1024/19*c - 1024/19, 1>,
//     <b^2 + 160/29*b*c - 320/29*b + 256/29*c^2 - 1024/29*c + 1024/29, 1>
// ]
Factorization(Parent(b)!(Evaluate(Coefficients(E)[5],b/(4*(c-2)))*(c-2)^6)); 
// [
//     <b^2 - 128/11*b*c + 256/11*b + 256/11*c^2 - 1024/11*c + 1024/11, 1>,
//     <b^2 + 160/29*b*c - 320/29*b + 256/29*c^2 - 1024/29*c + 1024/29, 2>
// ]
phi5facts := Factorization(DivisionPolynomial(E,5));
phi5facts[1];