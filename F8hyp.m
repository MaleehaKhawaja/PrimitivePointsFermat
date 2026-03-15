_<x>:=PolynomialRing(Rationals());
X:=HyperellipticCurve(-4*x^8+1);
J:=Jacobian(X);
lb,ub:=RankBounds(J);
assert lb eq 0 and ub eq 0; //J(Q) has rank 0
//Using the program of Reitsma and Muller, we find that #J(Q)=4 

D0:=Place(X![0,1,1]);
D1:=Place(X![0,-1,1]);

//Generator of J(Q)
assert IsPrincipal(D1-D0) eq false;
assert IsPrincipal(2*(D1-D0)) eq false;
assert IsPrincipal(4*(D1-D0)) eq true;

//------------------------------------------------//


for a in [1,2,3] do	
    D:=a*(D1-D0)+2*D0+2*D1;
    L,phi:=RiemannRochSpace(D);
    assert Dimension(L) eq 2;
end for;

//Case a = 0
F<u,v>:=FunctionField(X);
D:=0*(D1-D0)+2*D0+2*D1;
L,phi:=RiemannRochSpace(D);
assert Dimension(L) eq 3;
assert {phi(L.1), phi(L.2), phi(L.3)} eq {1/u^2,1/u,1}; //Quartic field generated is imprimitive
//------------------------------------------------//

//Case a = 2

Qu<u>:=FunctionField(Rationals());
Qux<x>:=PolynomialRing(Qu);

X:=HyperellipticCurve(-4*x^8+1);

D0:=Place(X![0,-1,1]);
D1:=Place(X![0,1,1]);

a:=2;

D:=a*(D1-D0)+2*D0+2*D1;
L,phi:=RiemannRochSpace(D);
assert Dimension(L) eq 2;

F<r,s>:=FunctionField(X);
phi(L.1);
phi(L.2);

// Plane model of C_8

Aff<x,y>:=AffineSpace(Rationals(), 2);

C1:=Curve(Aff, y^2 + 4*x^8 - 1);
C2:=Curve(Aff, 4*x^4 + y^2*x^4 - 2*y); //Plane model

C:=ProjectiveClosure(C1);
H:=ProjectiveClosure(C2);

F<a,b>:=FunctionField(C);

bool,phi:=IsIsomorphic(C,H);

assert bool eq true;

eqs:=DefiningEquations(phi);
eqsab:=[Evaluate(f,[a,b,1]): f in eqs];

A:=eqsab[1]/eqsab[3];
B:=eqsab[2]/eqsab[3];

assert 4*A^4+(B^2)*(A^4)-2*B eq 0;

//---------------------------------------------//
A<r,s>:=AffineSpace(Rationals(), 2);
//C2:=Curve(Aff, 4*x^4 + y^2*x^4 - 2*y); //Plane model

f:=4*r^4+r^4*s^2-2*s;

Discriminant(f, r);

_<t>:=PolynomialRing(Rationals());

H:=HyperellipticCurve(-2*t^3-8*t);
E:=EllipticCurve(H);

CremonaReference(E);
MordellWeilGroup(E);
