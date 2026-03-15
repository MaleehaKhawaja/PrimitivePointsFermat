//Verify the morphism from C_6 to the elliptic curve with Cremona label 27a3

_<x>:=PolynomialRing(Rationals());
C:=HyperellipticCurve(-4*x^6+1);
Aut:=Automorphisms(C);
w:=Aut[4];
G:=AutomorphismGroup(C, [w]);

E,phi:=CurveQuotient(G);

E;

assert Degree(phi) eq 2;

CremonaReference(E);

MordellWeilGroup(E);

eqs:=DefiningEquations(phi);

K<a,b>:=FunctionField(C);
eqsab:=[Evaluate(f,[a,b,1]): f in eqs];

A:=eqsab[1]/eqsab[3];
B:=eqsab[2]/eqsab[3];

assert B^2 + B - A^3 eq 0;

A;

B;