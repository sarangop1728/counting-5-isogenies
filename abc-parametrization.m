load "canring.m";

// We compute modular forms a,b, and c
// for the ring Gamma0(5) satisfying
// a^2 + b^2 = c^4.

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

// Check that a^2 + b^2 = c^4.
MonomialsOfWeightedDegree(R, 8);
Kernel(Transpose(Matrix([[Coefficient(g,n) : g in [a^2,b^2,c^4]]: n in [0..20]])));
a^2 + b^2 eq c^4;

// Equations for E4 and E6 in
// terms of a,b,c.

MonomialsOfWeightedDegree(R,4);
Kernel(Transpose(Matrix([[Coefficient(g,n) : g in [e4,a,b,c^2]]: n in [0..20]])));
e4 eq 246*a + 228*b + 250*c^2;

e6 := Eisenstein(6,q : Precision := 100);
MonomialsOfWeightedDegree(R,6);
Kernel(Transpose(Matrix([[Coefficient(g,n) : g in [e6,a*c,b*c,c^3]]: n in [0..20]])));
e6 eq -10008*a*c - 1044*b*c - 10000*c^3;
