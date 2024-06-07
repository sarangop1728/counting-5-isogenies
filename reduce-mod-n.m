P<a,b,c> := PolynomialRing(Integers(), 3);

// given a polynomial f in P, and an integer n, reduce the polynomial f modulo n.
function reduce(f, n)
    S := ChangeRing(P, Integers(n));
    h := hom< P -> S | a,b,c >;
    return h(f);
end function;

// example:
A := -6*(123*a + 114*b + 125*c^2);
B := 8*c*(2502*a + 261*b + 2500*c^2);
reduce(A,25);
