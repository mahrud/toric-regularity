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


