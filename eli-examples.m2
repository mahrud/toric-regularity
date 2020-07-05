-- Eli's examples
---attempt to do d=3 for fixed height a
-------------------------
restart
needsPackage "FourTiTwo"

--codim 2 ONLY right now (4 columns, 2 rows)
n = 5  --7 is the smallest integer where all four lists get populated
h = 3  --the fixed height


L = {}; ----making the list of all possible columns

for i to h-1 do (
    for j to h-1 do (
    	for k to h-1 do (
	    if i+j+k == h then (
	    m := matrix{{i},{j},{k}};
	    L = L|{m};
	    )
	)
    )
)

L


--Make lists we wish to populate:
CMeq = {};
noCMeq = {};
CMless = {};
noCMless = {};
--scan(length L',k->(
    	SL := subsets(L,2); --all subsets of columns we could append to our first column (1,0)^t, make this more general later....
	scan(SL,j->(
		K := matrix{{h, 0, 0}, {0,h,0}, {0,0,h}} ; --our first column
		scan(j,s-> K = K|s); --make a matrix K to consider
		S := QQ[x_1..x_(numcols K)];
		IA := toricMarkov(K,S); --make the toric ideal for the matrix K
		--Check that ZZA == ZZ^d, IA\subseteq (vars)^2, and IA nonzero:
		if ZZ^(numrows K) == image K then if IA != 0 then if (gens IA) % (ideal S_*)^2 == 0 then(
    		if regularity IA == (degree IA - codim IA + 1) then(
		    --we continue here when E-G holds with equality:
		    if pdim (S^1/IA) == codim IA then
		    CMeq = CMeq|{K} --CM case
		    else noCMeq = noCMeq|{K} --nonCM case
		    )
		else
		--we continue here when E-G does not hold with equality
		--(for curves, it will be strict inequality):
		    if pdim (S^1/IA) == codim IA then
		    CMless = CMless|{K} --CM case
		    else noCMless = noCMless|{K}; --nonCM case
		    )
    		))
--	))
--Here are the populated lists:
CMeq
CMless
noCMeq
noCMless

eqL = CMeq|noCMeq
noteqL = CMless|noCMless
