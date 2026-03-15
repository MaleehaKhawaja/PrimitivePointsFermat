//The code in this file supports the proof determining C_7(Q)

_<x>:=PolynomialRing(Rationals());
X:=HyperellipticCurve(-4*x^7+1);
J:=Jacobian(X);
lb,ub:=RankBounds(J);
assert lb eq 0 and ub eq 0; //J(Q) has rank 0
//Using the Magma program of Reitsma and Muller, we find that J(Q) = Z/7Z

assert Order(X![1,0,0]-X![0,-1,1]) eq 7;  //This is a generator of J(Q)

//We now determine all rational points on X

F<u,v>:=FunctionField(X);

D0:=Place(X![1,0,0]);
D1:=Place(X![0,-1,1]);

for a in [0..6] do	
    D:=a*(D0-D1) + D1;
    L,phi:=RiemannRochSpace(D);
    print "For a = ", a, "the RR dimension is", Dimension(L);
    if Dimension(L) eq 1
        then print phi(L.1);
    end if;
end for;

//Only rational points are [1,0,0], [0,-1,1] and [0,1,1]