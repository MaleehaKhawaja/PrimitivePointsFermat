_<x>:=PolynomialRing(Rationals());
C:=HyperellipticCurve(-4*x^8+1);
Aut:=Automorphisms(C);
w:=Aut[3];
G:=AutomorphismGroup(C, [w]);
E,phi:=CurveQuotient(G);

E; 
assert Degree(phi) eq 2;

CremonaReference(E); //64a4

eqs:=DefiningEquations(phi);
KC<a,b>:=FunctionField(C);
eqsab:=[Evaluate(f,[a,b,1]): f in eqs];

A:=eqsab[1]/eqsab[3];
B:=eqsab[2]/eqsab[3];

assert B^2 - A^3 - A eq 0;

A;

B;

