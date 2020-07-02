-- Mahrud's examples

-- Examples from June 30th
-- example from Nitsche paper

restart
alpha = 4
B = {{4,0},{0,4},{1,3},{3,1}}
R = QQ[a, b, c, d, Degrees => B]

k = 10
matrix table(k, k, (i, j) -> rank basis({k-i-1, j}, R))
select(alpha, i -> 0 == rank basis({i,alpha-i}, R))
apply(2+1, j -> select(alpha*j, i -> 0 == rank basis({i, alpha*j - i}, R)))
