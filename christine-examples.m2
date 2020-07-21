restart
needsPackage "FourTiTwo"

--Here is one way to compute a curve's arithmetic genus in M2:

K = matrix{{1,1,1,1},{0,1,3,4}}
S = QQ[x_1..x_(numcols K)]
IA = toricMarkov(K,S)
genus(integralClosure(Proj(S/IA)))

------
--Here is one way to compute a curve's geometric genus in M2:

K = matrix{{1,1,1,1},{0,1,3,4}}
S = QQ[x_1..x_(numcols K)]
IA = toricMarkov(K,S)
rank HH^1 sheaf(integralClosure(S/IA))

---------------
---------------

--Here is one way to compute a variety's arithmetic genus in M2:

K = matrix{{1,1,1,1},{0,1,3,4}}
S = QQ[x_1..x_(numcols K)]
IA = toricMarkov(K,S)
genus Proj(integralClosure(S/IA))

------
--Here is one way to compute a variety's geometric genus in M2:

K = matrix{{1,1,1,1},{0,1,3,4}}
S = QQ[x_1..x_(numcols K)]
IA = toricMarkov(K,S)
X = Proj(integralClosure(S/IA))
rank HH^0 ((cotangentSheaf(X))^(dim(X)))


----------------------------------------------------
----------------------------------------------------
----------------------------------------------------
--07/09/2020
--Attempts at understanding p.187 in Peeva--Sturmfels

restart
needsPackage "FourTiTwo"
needsPackage "Binomials"

A = matrix{{1,1,1,1,1},{0,1,0,1,0},{0,0,1,1,-2}}
S = QQ[x_1..x_(numcols A)]
IA = toricMarkov(A,S)
regularity IA
degree IA --equality holds in this case
B = gens kernel A

J = saturate(eliminate(IA+ideal(x_3-x_5),x_5),x_1*x_2*x_3*x_4)
regularity J
degree J
binomialAssociatedPrimes J
betti res J
betti res IA

J = saturate(eliminate(IA+ideal(x_4-x_5),x_5),x_1*x_2*x_3*x_4)
regularity J
degree J
binomialAssociatedPrimes J

B
B^{2}+B^{4}
B^{0}

B' = B^{0,1,2,3}
toricMarkov(transpose gens kernel transpose B',S)
B' = B^{0,1,2,5}
IA == toricMarkov(transpose gens kernel transpose B',S)

-----
restart
needsPackage "FourTiTwo"
needsPackage "Binomials"

A = transpose matrix{{1,0,0},{1,2,0},{1,3,0},{1,0,1},{1,1,1}}
S = QQ[x_1..x_(numcols A)]
IA = toricMarkov(A,S)
regularity IA
degree IA --equality holds in this case
B = gens kernel A

J = saturate(eliminate(IA+ideal(x_1-x_3),x_1),x_5*x_2*x_3*x_4) --I_{L'} on p.187 in PS
regularity J
degree J
binomialAssociatedPrimes J
betti res J
betti res IA

J = saturate(eliminate(IA+ideal(x_1-x_3),x_3),x_1*x_2*x_5*x_4)
regularity J
degree J
binomialAssociatedPrimes J
betti res J
betti res IA

J = saturate(eliminate(IA+ideal(x_3-x_5),x_5),x_1*x_2*x_3*x_4)
regularity J
degree J
binomialAssociatedPrimes J
betti res J
betti res IA

J = saturate(eliminate(IA+ideal(x_3-x_5),x_3),x_1*x_2*x_5*x_4)
regularity J
degree J
binomialAssociatedPrimes J
betti res J
betti res IA

-----
restart
needsPackage "FourTiTwo"
needsPackage "Binomials"

A = transpose matrix{{1,0,0},{1,1,0},{1,2,0},{1,1,1},{1,0,2}}
S = QQ[x_1..x_(numcols A)]
IA = toricMarkov(A,S)
regularity IA
degree IA --equality holds in this case
B = gens kernel A

J = saturate(eliminate(IA+ideal(x_1-x_3),x_1),x_5*x_2*x_3*x_4)
regularity J
degree J
binomialAssociatedPrimes J
betti res J
betti res IA

J = saturate(eliminate(IA+ideal(x_1-x_3),x_3),x_5*x_2*x_1*x_4)
regularity J
degree J
binomialAssociatedPrimes J
betti res J
betti res IA

J = saturate(eliminate(IA+ideal(x_1-x_5),x_1),x_5*x_2*x_3*x_4)
regularity J
degree J
binomialAssociatedPrimes J
betti res J
betti res IA


J = saturate(eliminate(IA+ideal(x_1-x_5),x_5),x_1*x_2*x_3*x_4)
regularity J
degree J
binomialAssociatedPrimes J
betti res J
betti res IA




-------------------------------
-------------------------------
--07/20/2020
--Potential PS counterexample? 
--No: there was a non-minimizing issue with M2
restart
needsPackage "FourTiTwo"
needsPackage "Binomials"

A = matrix{{1,1,1,1,1},{0,0,0,1,1},{0,5,7,4,5}}
entries transpose A
--gale_A = matrix{{0,1},{1,-2},{-1,1},{-2,-3},{2,3}}

S = QQ[x_1..x_(numcols A),Degrees=>entries transpose A]
--S = QQ[x_1..x_(numcols A)]
IA = toricMarkov(A,S)
regularity IA
degree IA --equality holds in this case
res IA
(res IA).dd
betti res IA
degrees (res IA)_3
--{9,2,44} is the "C" to use to get the syz quadrangle
--Here it is:
max degrees (res IA)_3

--Need: monomial x^a with Aa=C
a' = solve(A,transpose matrix {{9,2,44}})
B = gens kernel A
--Christine's HW: quickly find a in the fiber of C (done)
--hack to keep going: 
--a = a'+7*B_{0} --x_2^7x_4x_5
--Updated: 07/21/2020at 4:45pm: This will find a nonneg vector in the fiber of C
a = (latticePoints polyhedronFromHData(-map(ZZ^5,ZZ^5,1),transpose matrix {{0,0,0,0,0}}, A,transpose matrix {{9,2,44}}))#0

--Q: Is P_a the unit square?
needsPackage "Polyhedra"
--help polyhedronFromHData
--conv(u: (Bu\leq a))
Pa = polyhedronFromHData(B,a)
Latt = latticePoints Pa  
--Find the "middle" vectors and convert to a matrix: 
COC = matrix {toList (set Latt)-(set {sub(1/2*sub(sum Latt,QQ),ZZ),map(ZZ^2,ZZ^1,0)})}

--Double check: do we care if COC has determinant 1? 

Newb = B*COC
--Which two vectors are in the same quadrant? Add this code? 
 

--Just checking: 
latticePoints polyhedronFromHData(newB,a) --It's a square!


    
--variable reduction x_1=x_5 (first quadrant)
--regularity IA should be 7, while regularity after reduction is 6
J = ideal mingens saturate(eliminate(IA+ideal(x_1-x_5),x_5),x_1*x_2*x_3*x_4)
codim J
regularity J
degree J
binomialAssociatedPrimes J 
isPrime J
betti res J
betti res IA

J = ideal mingens saturate(eliminate(IA+ideal(x_1-x_3),x_3),x_1*x_2*x_5*x_4)
codim J
regularity J
degree J
binomialAssociatedPrimes J 
isPrime J
betti res J
betti res IA
