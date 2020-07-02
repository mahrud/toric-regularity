--Codim 2 curves:
--vijay's comment
restart
needsPackage "FourTiTwo"

--codim 2 ONLY right now (4 columns, 2 rows)
n = 9 --7 is the smallest integer where all four lists get populated
L = apply(n,i-> matrix{{1},{i}}) --list of all possible column vectors
L' = drop(L,1) --We will always choose (1,0)^t as the first column of every matrix.
--Make lists we wish to populate:
CMeq = {};
noCMeq = {};
CMless = {};
noCMless = {};
--scan(length L',k->( 
    	SL := subsets(L',3); --all subsets of columns we could append to our first column (1,0)^t
	scan(SL,j->(
		K := L#0; --our first column
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

apply(eqL,K->( 
    gens ker K
    ))



apply(eqL,K->( 
	S := QQ[x_1..x_(numcols K)];
	IA := toricMarkov(K,S);
	--{regularity IA, regularity (S^1/IA), degree IA, degree (S^1/IA)}
	mingens IA
    ))


--Check on Barrera's thesis
S = QQ[x_1..x_4]
IA = toricMarkov(eqL#0,S)
mingens IA
B = gens kernel eqL#0
B*transpose matrix{{1,0}}
B*transpose matrix{{1,1}}
B*transpose matrix{{3,2}}

eqL
noteqL


apply(eqL,K->( 
	S := QQ[x_1..x_(numcols K)];
	IA := toricMarkov(K,S);
	peek res IA
    ))
apply(noteqL,K->( 
	S := QQ[x_1..x_(numcols K)];
	IA := toricMarkov(K,S);
	betti res IA
    ))


apply(eqL,K->( 
	S := QQ[x_1..x_(numcols K)];
	IA := toricMarkov(K,S);
	genus(Proj(S/IA))
    ))
apply(noteqL,K->( 
	S := QQ[x_1..x_(numcols K)];
	IA := toricMarkov(K,S);
	genus(Proj(S/IA))
    ))

------------
-- 1. Curves
-- Gruson--Lazarsfeld--Peskine: proved equality for curves! 
-- Use their conditions to translate to toric case! 
-- Goal: combinatorial description for A that determines 
-- when equality holds. 
-- Next steps: read the statements of GLP paper, understand 
-- their logical implications. Then determine which statement
-- should be translated into combinatorics of A. 

-- 2. Start an overleaf document to begin tabulating cases 
-- we understand. 

-- 3. d=3, n=5: Understand this case next! 
-- (Start with M2 loop. Find nice way to search!) 









-------------------------
-------------------------
--does the genus distinction generalize to all monomial curves? --No!

restart
needsPackage "FourTiTwo"

--codim 2 ONLY right now (4 columns, 2 rows)
n=9 --7 is the smallest integer where all four lists get populated
L = apply(n,i-> matrix{{1},{i}}) --list of all possible column vectors
L' = drop(L,1) --We will always choose (1,0)^t as the first column of every matrix.
--Make lists we wish to populate:
CMeq = {};
noCMeq = {};
CMless = {};
noCMless = {};
--scan(length L',k->( 
    	SL := subsets(L',4); --all subsets of columns we could append to our first column (1,0)^t
	scan(SL,j->(
		K := L#0; --our first column
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

apply(eqL,K->( 
	S := QQ[x_1..x_(numcols K)];
	IA := toricMarkov(K,S);
	genus(Proj(S/IA))
    ))
apply(noteqL,K->( 
	S := QQ[x_1..x_(numcols K)];
	IA := toricMarkov(K,S);
	genus(Proj(S/IA))
    ))

noteqL#24
