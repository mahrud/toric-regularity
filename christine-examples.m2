restart
needsPackage "FourTiTwo"

--Here is one way to compute a curve's arithmetic genus in M2:

K = matrix{{1,1,1,1},{0,1,3,4}}
S = QQ[x_1..x_(numcols K)]
IA = toricMarkov(K,S)
genus(Proj(S/IA))

------
--Here is one way to compute a curve's geometric genus in M2:

K = matrix{{1,1,1,1},{0,1,3,4}}
S = QQ[x_1..x_(numcols K)]
IA = toricMarkov(K,S)
rank HH^1 sheaf(S/IA)

---------------
---------------

--Here is one way to compute a variety's arithmetic genus in M2:

K = matrix{{1,1,1,1},{0,1,3,4}}
S = QQ[x_1..x_(numcols K)]
IA = toricMarkov(K,S)
genus Proj(S/IA)

------
--Here is one way to compute a variety's geometric genus in M2:

K = matrix{{1,1,1,1},{0,1,3,4}}
S = QQ[x_1..x_(numcols K)]
IA = toricMarkov(K,S)
X = Proj(S/IA)
rank HH^0 ((cotangentSheaf(X))^(dim(X)))
