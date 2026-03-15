_<x>:=PolynomialRing(Rationals());
X:=HyperellipticCurve(-4*x^7+1);
J:=Jacobian(X);

D0:=Place(X![1,0,0]);
D1:=Place(X![0,-1,1]);

for a in [1..6] do	
    D:=a*(D1-D0)+4*D0;
    L,phi:=RiemannRochSpace(D);
    assert Dimension(L) eq 2;
end for;

//Case a = 0
F<u,v>:=FunctionField(X);
D:=4*D0;
L,phi:=RiemannRochSpace(D);
assert Dimension(L) eq 3;

assert {phi(L.1), phi(L.2), phi(L.3)} eq {u^2,u,1}; //Quartic field generated is imprimitive

//-----------------------------------------------------------------------------------//

//Verifying plane model of C_7 for a = 3

Qu<u>:=FunctionField(Rationals());
Qux<x>:=PolynomialRing(Qu);

X:=HyperellipticCurve(-4*x^7+1);

D0:=Place(X![1,0,0]);
D1:=Place(X![0,-1,1]);

a:=3;

D:=a*(D1-D0)+4*D0;
L,phi:=RiemannRochSpace(D);
assert Dimension(L) eq 2;

F<m,n>:=FunctionField(X);
phi(L.1);
phi(L.2);

// Plane model of C_7

Aff<x,y>:=AffineSpace(Rationals(), 2);

C1:=Curve(Aff, y^2 + 4*x^7 - 1);
C2:=Curve(Aff, 4*x^4 + y^2*x^3 + 2*y); //Plane model

C:=ProjectiveClosure(C1);
H:=ProjectiveClosure(C2);

F<a,b>:=FunctionField(C);

bool,phi:=IsIsomorphic(C,H);

assert bool eq true;

eqs:=DefiningEquations(phi);
eqsab:=[Evaluate(f,[a,b,1]): f in eqs];

A:=eqsab[1]/eqsab[3];
B:=eqsab[2]/eqsab[3];

assert 4*A^4+(B^2)*(A^3)+2*B eq 0;


//----------------------------------------------------------------//

//Verifying plane model of C_7 for a = 4

Qu<u>:=FunctionField(Rationals());
Qux<x>:=PolynomialRing(Qu);

X:=HyperellipticCurve(-4*x^7+1);

D0:=Place(X![1,0,0]);
D1:=Place(X![0,-1,1]);

a:=4;

D:=a*(D1-D0)+4*D0;
L,phi:=RiemannRochSpace(D);
assert Dimension(L) eq 2;

F<m,n>:=FunctionField(X);
phi(L.1);
phi(L.2);

Aff<x,y>:=AffineSpace(Rationals(), 2);

C1:=Curve(Aff, y^2 + 4*x^7 - 1);
C2:=Curve(Aff, y^2*x^4 + 4*x^3 - 2*y); //Plane model

C:=ProjectiveClosure(C1);
H:=ProjectiveClosure(C2);

F<a,b>:=FunctionField(C);

bool,phi:=IsIsomorphic(C,H);

assert bool eq true;

eqs:=DefiningEquations(phi);
eqsab:=[Evaluate(f,[a,b,1]): f in eqs];

A:=eqsab[1]/eqsab[3];
B:=eqsab[2]/eqsab[3];

assert B^2*A^4+4*(A^3)-2*B eq 0;

//-----------------------------------------------------------------//

_<x>:=PolynomialRing(Rationals());
H:=HyperellipticCurve(2*x^7-27);

//Want to show H(Q) = {pt at infinity}

J:=Jacobian(H);
ub,lb:=RankBounds(J);
assert ub eq 0 and lb eq 0;

assert BadPrimes(J) eq [2,3,7];

bdp:=[2,3,7];

g0:=#BaseChange(J, GF(5));
count:=1;

for p in [11..50]
    do if IsPrime(p) 
        then count := count + 1;
            Jp:=#BaseChange(J, GF(p));
            g0:=Gcd(g0, Jp);
            if g0 eq 1
                then print "The MW group of the Jacobian is trivial";
                break;
            end if;
        end if;
end for;

print count;


