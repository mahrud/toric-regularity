
--Here is one way to compute genus in M2:

K = {{1,1,1,1},{0,1,3,4}}
S = QQ[x_1..x_(numcols K)]
IA = toricMarkov(K,S)
genus(Proj(S/IA))
