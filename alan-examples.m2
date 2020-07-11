--This program computes a list (namely, eqL) of certain 3xn matrices A such that I_A satisfies the conditions of EG and gives equality. We standardize our matrix A so that it has all nonnegative entries, its first row is all 1's, and all other rows start with 0's (it's not clear if this is the most natural way to standardize, but at least it reduces computation load). In the below computation we require all other entries to be nonnegative integers less than N, and give results up to a permutation of the non-first columns of A (note we assume A has distinct columns, as otherwise I_A does not satisfy the conditions). I would replace 3xn with the more general dxn but I'm too lazy to figure out how.

restart
needsPackage "FourTiTwo"

n = 5 --number of columns of A (we care about 3x5 the most right now)

N = 4 --put like 4 or something whatever floats your boat (but not too big, like 7)


--Make lists we wish to populate:
eqL = {};
noteqL = {}; --strictly speaking, this should be "less than" and not "not equal"
bad = {}; --A such that I_A doesn't satisfy the conditions of EG

CMeq = {};
noCMeq = {};
CMless = {};
noCMless = {};
notsurjective = {};
surjectivebutnotinsquareideal = {};
yeet = {};

--The below function takes a given matrix A and sorts it into the correct lists (note that this "function" is bad because it only works if the lists eqL, noteqL, etc. already exist, but whatever).
test = A->(
    if ZZ^(numrows A) != image A then(notsurjective = notsurjective|{A}; bad = bad|{A})--Check that ZZA == ZZ^d
    else(
	S := QQ[x_1..x_(numcols A)];
	IA := toricMarkov(A,S); --make the toric ideal for the matrix A
	--Check that IA\subseteq (vars)^2 (I removed checking that IA is zero): 
	if (gens IA) % (ideal S_*)^2 != 0 then (surjectivebutnotinsquareideal = surjectivebutnotinsquareideal|{A}; bad = bad|{A})
	else(
	    if regularity IA == (degree IA - codim IA + 1) then(
	    	    --we continue here when E-G holds with equality:
		    eqL = eqL|{A};
	    	    if pdim (S^1/IA) == codim IA then CMeq = CMeq|{A} --CM case
	    	    else noCMeq = noCMeq|{A} --nonCM case
		    )
    	    else(
		if regularity IA > (degree IA - codim IA + 1) then(yeet = yeet|{A}) --yeet case
		else(
		    noteqL = noteqL|{A};
		    if pdim (S^1/IA) == codim IA then(CMless = CMless|{A}) --CM case
		    else noCMless = noCMless|{A} --nonCM case
		    )
		)
	    )
	)
    )

L = {}; --this will soon become our list of non-first columns of A
scan(N,i-> (
	scan(N,j-> (L = L|{matrix{{1},{i},{j}}}))
	)
    )
L = drop(L,1); --this is now our list of non-first columns of A

SL := subsets(L,n-1); --this is a list of all possible combinations of non-first columns of A

scan(SL, j->(
	A := matrix{{1},{0},{0}}; -- this is the first column of A
	scan(j,s-> A = A|s); --this builds the rest of A, column by column
	test(A)
	)
    )

--Here are the populated lists:
eqL
CMeq
noCMeq

noteqL;
CMless;
noCMless;

bad;
notsurjective;
surjectivebutnotinsquareideal;

yeet --if this is nonempty life will get real interesting real fast