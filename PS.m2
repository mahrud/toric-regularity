--The below is similar to the stuff in alan-examples, but in the codim-2 non-CM case, will print out data associated with the reduction on pg. 187 of PS.

restart
needsPackage "FourTiTwo"
needsPackage "Polyhedra"



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
yeet = {}; --will be empty for codim 2...hopefully




gale = A->( --ONLY WORKS FOR NON-CM CODIM 2
    assert ((numcols A)-2 == numrows A);
    S := QQ[x_1..x_(numcols A),Degrees=>entries transpose A];
    IA := toricMarkov(A,S);
    
    --Here is the C we use to get the syzygy quadrangle:
    C := max degrees (res IA)_3;
    
    --Construct a with nonnegative entries such that Aa=C:
    zeros := {};
    scan(numcols A, i->(zeros = zeros|{0}));
    a := (latticePoints polyhedronFromHData(-map(ZZ^(numcols A),ZZ^(numcols A),1),transpose matrix {zeros}, A,transpose matrix {C}))#0;
    
    B := gens ker A;
    Pa := polyhedronFromHData(B,a);
    Latt := latticePoints Pa;
    
    --Find the "middle" vectors and convert to a matrix:
    COC := matrix {Latt-(set {sub(1/2*sub(sum Latt,QQ),ZZ),map(ZZ^2,ZZ^1,0)})};
    
    B*COC
    )





reduction = A->( --ONLY WORKS FOR NON-CM CODIM 2. ALSO SORRY FOR BAD CODE
    B := gale(A);
    S := QQ[x_1..x_(numcols A)];
    IA := toricMarkov(A,S);
    
    --Computes all possible reductions:
    R := QQ[y_1..y_4];
    images := {{}};
    scan(numrows B, i->(
	    imagesnew := {};
	    if B_(i,1)>0 then(
		if B_(i,0)>0 then(
		    scan(images,j->(imagesnew = imagesnew|{j|{R_0}}))
		    )
		else(
		    if B_(i,0)==0 then(
			scan(images,j->(imagesnew = imagesnew|{j|{R_0}}));
			scan(images,j->(imagesnew = imagesnew|{j|{R_1}}))
			)
		    else(
			scan(images,j->(imagesnew = imagesnew|{j|{R_1}}))
			)
		    )
		)
	    else(
		if B_(i,1)==0 then(
		    if B_(i,0)>0 then(
			scan(images,j->(imagesnew = imagesnew|{j|{R_0}}));
			scan(images,j->(imagesnew = imagesnew|{j|{R_3}}))
			)
		    else(
			if B_(i,0)==0 then(
			    scan(images,j->(imagesnew = imagesnew|{j|{R_0}}));
			    scan(images,j->(imagesnew = imagesnew|{j|{R_1}}));
			    scan(images,j->(imagesnew = imagesnew|{j|{R_2}}));
			    scan(images,j->(imagesnew = imagesnew|{j|{R_3}}))
			    )
			else(
			    scan(images,j->(imagesnew = imagesnew|{j|{R_1}}));
			    scan(images,j->(imagesnew = imagesnew|{j|{R_2}}))
			    )
			)
		    )
		else(
		    if B_(i,0)>0 then(
			scan(images,j->(imagesnew = imagesnew|{j|{R_3}}))
			)
		    else(
			if B_(i,0)==0 then(
			    scan(images,j->(imagesnew = imagesnew|{j|{R_2}}));
			    scan(images,j->(imagesnew = imagesnew|{j|{R_3}}))
			    )
			else(
			    scan(images,j->(imagesnew = imagesnew|{j|{R_2}}))
			    )
			)
		    )
		);
	    images = imagesnew
	    )
	);
    
    apply(images,i->saturate((map(R, S, i))(IA), ideal(R_0*R_1*R_2*R_3)))
    )




--Examples to try reduction on:
A = matrix{{1, 1, 1, 1, 1}, {0, 0, 0, 1, 1}, {0, 5, 7, 4, 5}} -- (3x5, non-CM, equality; this is our "favorite" example)
A1 = matrix{{1, 1, 1, 1, 1, 1}, {0, 0, 0, 0, 1, 1}, {0, 0, 1, 1, 0, 0}, {0, 1, 1, 2, 0, 2}} -- (4x6, non-CM, equality)
A2 = matrix{{1, 1, 1, 1, 1, 1}, {0, 0, 0, 1, 1, 1}, {0, 1, 1, 0, 1, 1}, {0, 1, 2, 2, 0, 1}} -- (4x6, non-CM, less)






--testPS takes a given matrix A and sorts it into the correct list, but in the codimension-2 non-CM case, it also outputs the Gale diagram from PS pg. 187 (more precisely, the matrix B), along with the regularity, degree, and Betti table corresponding to A and each of its reductions. Note the Gale diagram is *not necessarily* uniquely determined from A, so other reductions may be possible.



PSeq = {};
PSless = {};

testPS = A->(
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
	    	    else(
			noCMeq = noCMeq|{A}; --nonCM case
			if codim IA == 2 then( --PS reduction for nonCM and codim 2
			    reductions := {};
			    scan(reduction(A), IL'->(reductions = reductions|{{regularity IL', degree IL', betti res IL'}}));
			    PSeq = PSeq|{{A, gale(A), {regularity IA, degree IA, betti res IA}, reductions}}
			    )
			)
		    )
    	    else(
		if regularity IA > (degree IA - codim IA + 1) then(yeet = yeet|{A}) --yeet case
		else(
		    noteqL = noteqL|{A};
		    if pdim (S^1/IA) == codim IA then(CMless = CMless|{A}) --CM case
		    else(
			noCMless = noCMless|{A}; --nonCM case
			if codim IA == 2 then( --PS reduction for nonCM and codim 2
			    reductions := {};
			    scan(reduction(A), IL'->(reductions = reductions|{{regularity IL', degree IL', betti res IL'}}));
			    PSless = PSless|{{A, gale(A), {regularity IA, degree IA, betti res IA}, reductions}}
			    )
			)
		    )
		)
	    )
	)
    )






--Basic loops:



--3x5: feel free to change the {0,1,2,3} to other lists that start with 0

L = {}; --this will soon become our list of non-first columns of A
scan({0,1,2,3},i-> (
	scan({0,1,2,3},j-> (L = L|{matrix{{1},{i},{j}}}))
	)
    )
L = drop(L,1); --this is now our list of non-first columns of A

SL := subsets(L,4); --this is a list of all possible combinations of non-first columns of A

scan(SL, j->(
	A := matrix{{1},{0},{0}}; -- this is the first column of A
	scan(j,s-> A = A|s); --this builds the rest of A, column by column
	testPS(A)
	)
    )




--4x6: again, feel free to change the lists {0,1},{0,1},{0,1,2} to other lists that start with 0

L = {}; --this will soon become our list of non-first columns of A
scan({0,1},i-> (
	scan({0,1},j-> (
		scan({0,1,2},k-> (
			L = L|{matrix{{1},{i},{j},{k}}}
			)
		    )
		)
	    )
	)
    )
L = drop(L,1); --this is now our list of non-first columns of A

SL := subsets(L,5); --this is a list of all possible combinations of non-first columns of A

scan(SL, j->(
	A := matrix{{1},{0},{0},{0}}; -- this is the first column of A
	scan(j,s-> A = A|s); --this builds the rest of A, column by column
	testPS(A)
	)
    )