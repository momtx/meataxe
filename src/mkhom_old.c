/* ============================= C MeatAxe ==================================
   File:        $Id: mkhom_old.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Calculate homomorphisms between modules.
   --------------------------------------------------------------------------
   Written by Magdolna Szoke. 
   Revised by Michael Ringe.
   (C) Copyright 1999 Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"
#include <string.h>



MTX_DEFINE_FILE_INFO

static MtxApplicationInfo_t AppInfo = { 
"mkhom", "Calculate homomorphisms",
"SYNTAX\n"
"    mkhom [-ts] [-r <Side>] [-b <Mode>] [-H <Dim>] <M> <N> <Hom>\n"
"\n"
"ARGUMENTS\n"
"    <M> ..................... First representation\n"
"    <N> ..................... Second representation\n"
"    <Hom> ................... Homomorhisms from <M> to <N>\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -t ...................... Calculate standard generators for <M>\n"
"    -s ...................... When <M>=<N>, give endomorphisms in standard basis\n"
"    -r <Side>................ When <M>=<N>, find a generating set of End(M), and\n"
"                              calculate the left (<Side>=1) or right (<Side>=2)\n"
"                              regular representation.\n"
"    -b <Mode>................ Save memory, <Mode>=0..2.\n"
"    -H <Dim> ................ If the radical is given, <Dim> is the dimension of\n"
"\n"
"FILES\n"
"    <M>.{1,2...} ............ I  Generators in representation <M>.\n"
"    <N>.{1,2...} ............ I  Generators in representation <N>.\n"
"    <M>.cfinfo .............. I  Constituent info file for <M>.\n"
"    <N>.cfinfo .............. I  Constituent info file for <N>.\n"
"    <M>.rad ................. I  Generators for the head of <M> (with -H).\n"
"    <M><Cf>.k ............... I  Uncondense matrix, produced by PWKOND.\n"
"    <M>.std.................. O  The standard basis for <M>.\n"
"    <Hom>.{1,2,...} ......... O  A k-basis of Hom(<M>,<N>).\n"
"    <M>.std.{1,2,...} ....... O  Generators in standard basis (with -t).\n"
};

static MtxApplication_t *App = NULL;
static int standard = 0;		/* -t: Make standard generators */
static int hominstd = 0;		/* -s: Endomorphisms in std. basis */
static int reg = 0;			/* -r: Regular representation */
static char side = '?';			/* 'l' or 'r' (with -r) */
static int big = 0;			/* Argument of -b */
static int hd = 0;			/* Head dimension (from -H) */
static const char *TempDirName = NULL;	/* Temporary directory (with -b) */
static const char *MName = NULL;	/* First argument */
static const char *NName = NULL;	/* Second argument */
static const char *HomName = NULL;	/* Third argument */
static int comp = 0;			/* strcmp(MName,NName) */
static Lat_Info MInfo;			/* Constituent info */
static MatRep_t *MRep = NULL;		/* First representation */
static MatRep_t *NRep = NULL;		/* Second representation */
static int dim;				/* Dimension of <M> */

static PTR basis = NULL;
static PTR space = NULL;
static int partdim = 0;
static int *piv = NULL;
static long *op = NULL;
static PTR *stdgen = NULL;
static long **stdtab = NULL;
static long *tab = NULL;
static Matrix_t **oldstdbas = NULL;
static Matrix_t *rad = NULL;




/* TODO: ersetzen durch FfCleanRowAndRepeat */

static void myzcleanrow(PTR row, PTR matrix, PTR matrix2, int nor, int *piv)

{
    long i;
    PTR x, y, row2;

    row2 = FfGetPtr(matrix2,nor);
    for (i = 0, x = matrix, y = matrix2; i < nor; 
         ++i, FfStepPtr(&x), FfStepPtr(&y))
    {
        FEL f = FfExtract(row, piv[i]);
        if (f == FF_ZERO) 
	    continue;
        f = FfNeg(FfDiv(f, FfExtract(x, piv[i])));
        FfAddMulRow(row, x, f);
	FfAddMulRow(row2, y, f);
    }
}


/* ------------------------------------------------------------------
   zgensbasis() - Spin up canonically (standard basis).

   <seed> ist ein Zeiger auf Seedvektoren, Der seedcount - 1-ste
   wird benutzt. <gen> sind die Erzeuger (gen[0]...gen[ngen-1]).

   <space> und <basis> muessen vom Aufrufer allokiert werden. Beide
   muessen gross genug fuer eine quadratrische Matrix sein. In <basis>
   wird die Standardbasis abgelegt.

   <piv> ist die Pivot-Tabelle, die ebenfalls vom Aufrufer allokiert
   werden muss (mindestans dim + 1 Elemente!).

   In <op_table> wird die Definition der Standardbasisvektoren abgelegt.
   (2 * (dim+1) long integers). Falls diese Information nicht
   benoetigt wird, kann op_table=NULL gesetzt werden.

   das Stueck, das der <seedcount>-ste Vektor ueber den schon existierenden
   <partdim> dimensionalen Modul erzeugt (der von den ersten seedcount - 1
   Vektoren aufgespannt ist) wird berechnet.

   <stdgen> sind die Standard erzeuger (stueckweise aufgebaut),
   in <stdtab> steht, welche Operationen zu ueberpruefen sind um zu entscheiden,
   ob ein gewisser Modul mit aehnlicher Basis zu dem gegebenen Modul isomorph
   ist.
   <stdgen> muss fuer ngen PTR allokiert sein, sowie stdgen[i] fuer 0 Zeilen;
   <stdtab> muss fuer ngen * 1 long integers mit stdtab[i][0] == 0.

   Koerperordnung und Zeilenlaenge muessen vor dem Aufruf gesetzt
   werden.

   Return: neue Dimension.
   ------------------------------------------------------------------ */
#define OPVEC(i) op_table[2*(i)]
#define OPGEN(i) op_table[2*(i)+1]


int zgensbasis(PTR seed, int seedcount, int ngen, Matrix_t *gen[], PTR space, 
    int *piv_table, PTR basis, int partdim, long *op_table, PTR stdgen[],
    long *std_tab[])

{
    static Matrix_t *transf = NULL;
    static long gencount = 1;
    long i, j, k;
    PTR xi, yi, xk, yk, temp, row, transfptr;
    FEL f;
    int igen;

    /* Initialize
       ---------- */
    if (transf == NULL)	/* identity matrix + one zero-row for the operation */
    {
	transf = MatAlloc(FfOrder, FfNoc + 1, FfNoc);
	for (i = 0, row = transf->Data; i < FfNoc; ++i)
	{
	    FfInsert(row, i, FF_ONE);
	    FfStepPtr(&row);
	}
    }
    i = 1; 
    j = partdim + 1; 
    xi = FfGetPtr(space,partdim);
    yi = FfGetPtr(basis,partdim);
    k = partdim + 1; 
    xk = FfGetPtr(space,partdim);
    yk = FfGetPtr(basis,partdim);
    igen = 0;
    seed = FfGetPtr(seed,seedcount-1);

    /* Main loop
       --------- */
    FfCopyRow(yk,seed);
    FfCopyRow(xk,seed);
    if (op_table != NULL)
    {
        OPVEC(k) = gencount;
	OPGEN(k) = 0;
    }
    myzcleanrow(xk, space, transf->Data, partdim, piv_table);
    if ((piv_table[partdim] = FfFindPivot(xk, &f)) < 0)
    {
	transfptr = MatGetPtr(transf,partdim);
	FfMulRow(transfptr,FF_ZERO);
	if (partdim < FfNoc) 
	    FfInsert(transfptr, partdim, FF_ONE);
	return partdim;
    } 
    gencount++;
    k++;
    partdim++;
    FfStepPtr(&xk);
    FfStepPtr(&yk);

    while (xi < xk)
    {
	FfMapRow(yi, gen[igen]->Data, FfNoc, yk);
	FfCopyRow(xk,yk);

	/* Clean and check if we got a new vector
	   -------------------------------------- */
	myzcleanrow(xk, space, transf->Data, partdim, piv_table);
	if ((piv_table[partdim] = FfFindPivot(xk, &f)) >= 0)
	{
	    if (op_table != NULL)
	    {
		OPVEC(k) = j;
		OPGEN(k) = igen+1;
	    }
	    ++k;
	    partdim++;
	    FfStepPtr(&xk);
	    FfStepPtr(&yk);
	}
	else
	{
	    std_tab[igen] = NREALLOC(std_tab[igen],long,++std_tab[igen][0] + 1);
	    std_tab[igen][std_tab[igen][0]] = i;
	    temp = FfAlloc(std_tab[igen][0]);
	    memcpy(temp,stdgen[igen],FfCurrentRowSize * (std_tab[igen][0] - 1));
	    row = FfGetPtr(temp,std_tab[igen][0] - 1);
	    transfptr = FfGetPtr(transf->Data,partdim);
	    FfCopyRow(row,transfptr);
	    if (partdim < FfNoc) 
		FfInsert(row, partdim, FF_ZERO);
	    FfMulRow(row, FfNeg(FF_ONE));
	    FfMulRow(transfptr, FF_ZERO);
	    if (partdim < FfNoc) 
		FfInsert(transfptr, partdim, FF_ONE);
	    SysFree(stdgen[igen]);
	    stdgen[igen] = temp;
	}

	if (++igen >= ngen)	/* All the generators have been used */
	{
	    igen = 0;
	    ++i;
	    ++j;
	    FfStepPtr(&xi);
	    FfStepPtr(&yi);
	}

    }
    return partdim;
}


long independent(Matrix_t *bas[], Matrix_t *mat, int dim, int **piv_table, 
    const long *table, PTR dep)

{
    int i, j;
    FEL f;
    PTR matptr, basptr;

    MESSAGE(1,("independent: dim=%d\n",dim));
    FfSetNoc(mat->Noc);
    for (i = 0; i < dim; i++)
    {
	if (bas[i] == NULL)
	    continue;
	basptr = MatGetPtr(bas[i],piv_table[i][0]);
	matptr = MatGetPtr(mat,piv_table[i][0]);
	f = FfExtract(matptr, piv_table[i][1]);
	f = FfDiv(f, FfExtract(basptr, piv_table[i][1]));
	if (dep != NULL)
#ifdef PARANOID
	    FfSetNoc(dim);
#endif
	    FfInsert(dep,i,f);
	MatAddMul(mat,bas[i],FfNeg(f));
    }
    piv_table[dim][0] = -1;
    if (table == NULL)
    {
	matptr = mat->Data;
	for (j = 0; j < mat->Nor && piv_table[dim][0] < 0; j++)
	{
	    if ((piv_table[dim][1] = FfFindPivot(matptr, &f)) >= 0)
		piv_table[dim][0] = j;
	    FfStepPtr(&matptr);
	}
    }
    else
    {
	int row = 0;
	matptr = mat->Data;
	for (j = 1; j <= table[0] && piv_table[dim][0] < 0; j++)
	{
	    if ((piv_table[dim][1] = FfFindPivot(matptr, &f)) >= 0)
		piv_table[dim][0] = row;
	    matptr = FfGetPtr(matptr,table[j]);
	    row += table[j];
	}
    }

#ifdef PARANOID
    if (dim >= FfNoc) FfSetNoc(dim+1);
#endif
    if (piv_table[dim][0] >= 0 && dep != NULL)
	FfInsert(dep, dim, FF_ONE);
    MESSAGE(2,("independent(): result=%d\n",piv_table[dim][0] >= 0));
    return piv_table[dim][0] >= 0;
}








Matrix_t *bigform(Matrix_t *mat, Matrix_t **gens, long *op_table, int dim)

{
    long ind, max;
    Matrix_t *big;
    PTR bigptr, matptr, ptr;

    matptr = mat->Data;
    big = MatAlloc(mat->Field, dim, mat->Noc);
    bigptr = big->Data;
    max = 2 * dim;
    for (ind = 2; ind <= max; ind += 2)
    {
	if (op_table[ind + 1] == 0)
	{
	    FfCopyRow(bigptr, matptr);
	    FfStepPtr(&matptr);
	}
	else
	{
	    ptr = MatGetPtr(big,op_table[ind] - 1);
	    FfMapRow(ptr, gens[op_table[ind + 1] - 1]->Data, gens[0]->Nor, bigptr);
	}
	FfStepPtr(&bigptr);
    }
    return big;
}













Matrix_t **ringgens(Matrix_t *basis[], long n, long *table, Matrix_t *regrep[], 
    char side, int big, Matrix_t **stdbas, long *op_table, Matrix_t **Ngen)

{
    int max = 0, next = 0, i, j, dim = 0, 
	**piv_table, **bpiv, **baspiv, d, g, *genind, a, b, c;
    Matrix_t **gens, *mat;
    PTR *regptr;
    FEL coeff;

    if (side != 'l' && side != 'r')
    {
	MTX_ERROR1("Invalid side='%c'",side);
	return NULL;
    }
    if (   (piv_table = NALLOC(int *,n + 1)) == NULL
        || (baspiv = NALLOC(int *,n + 1)) == NULL)
	return NULL;

    for (i = 0; i <= n; i++)
    {
	if (   (piv_table[i] = NALLOC(int,2)) == NULL
            || (baspiv[i] = NALLOC(int,2)) == NULL)
	    return NULL;
    }
    if (   (genind = NALLOC(int,n)) == NULL
        || (gens = NALLOC(Matrix_t *,n + 1)) == NULL
        || (regptr = NALLOC(PTR,n)) == NULL
        || (bpiv = NALLOC(int *,2)) == NULL
        || (bpiv[0] = NALLOC(int,2)) == NULL
        || (bpiv[1] = NALLOC(int,2)) == NULL)
	return NULL;

/* -----------------------------
   makes a basis for the algebra
   ----------------------------- */

    d = basis[0]->Noc;
    g = basis[0]->Nor;
    while (dim < n)
    {
    	MESSAGE(1,("ringgens(): dim=%d\n",dim));
	if ((stdbas[dim] = MatAlloc(FfOrder, g, d)) == NULL)
	    return NULL;

/* chosing a random element of the algebra
   --------------------------------------- */
	for (i = 0; i < n; i++)
	{
	    coeff = FfFromInt(MtxRandomInt(FfOrder));
	    if (basis[i] != NULL)
		MatAddMul(stdbas[dim], basis[i], coeff);
	}

/* testing if it is independent from the others
   -------------------------------------------- */
	if (!independent(stdbas, stdbas[dim], dim, piv_table, table, NULL))
	{
	    MatFree(stdbas[dim]);
	    continue;
	}
	genind[max] = dim;
	switch (big)
	{
	    case 0:
		gens[max] = stdbas[dim];
		break;
	    case 1:
	    case 2:
	    case 3:
		gens[max] = bigform(stdbas[dim], Ngen, op_table, d);
		break;
	}
	if (big)
	{
	    bpiv[0][0] = piv_table[dim][0];
	    bpiv[0][1] = piv_table[dim][1];
	    c = 0;
	    for (i = 0; i < n; i++)
	    {
		if (basis[i] == NULL)
		    continue;
		a = independent(&(stdbas[dim]), basis[i], 1, bpiv, NULL, NULL);
		if (c) continue;
		b = independent(basis, basis[i], i, baspiv, NULL, NULL);
		if (!a || !b)
		{
		    MatFree(basis[i]);
		    basis[i] = NULL;
		    c = 1;
		}
	    }
	}
	dim++;
    	MESSAGE(1,("ringgens(): new element, dim=%d\n",dim));
	if ((regrep[max] = MatAlloc(FfOrder, n, n)) == NULL)
	    return NULL;
	regptr[max] = regrep[max]->Data;
	for (i = 0; i < genind[max]; i++)
	{
	    if (side == 'r')
	    {
		stdbas[dim] = MatDup(stdbas[i]);
		mat = stdbas[dim];
		MatMul(mat, gens[max]);
	    }
	    else 	/* side == 'l' -- multiple from left */
	    {
		stdbas[dim] = MatDup(stdbas[genind[max]]);
		if (big)
		    mat = bigform(stdbas[i], Ngen, op_table, d);
		else
		    mat = MatDup(stdbas[i]);
		MatMul(stdbas[dim], mat);
		MatFree(mat);
		mat = stdbas[dim];
	    }
	    if (independent(stdbas, mat, dim, piv_table, table, regptr[max]))
	    {
		if (big)
		{
		    bpiv[0][0] = piv_table[dim][0];
		    bpiv[0][1] = piv_table[dim][1];
		    c = 0;
		    for (j = 0; j < n; j++)
		    {
			if (basis[j] == NULL)
			    continue;
			a = independent(&(stdbas[dim]), basis[j], 1, bpiv, NULL, NULL);
			if (c) continue;
			b = independent(basis, basis[j], j, baspiv, NULL, NULL);
			if (!a || !b)
			{
			    MatFree(basis[j]);
			    basis[j] = NULL;
			    c = 1;
			}
		    }
		}
		dim++;
    		MESSAGE(1,("ringgens(): new element2, dim=%d\n",dim));
	    }
	    else
		MatFree(stdbas[dim]);

	    FfSetNoc(n);
	    FfStepPtr(&(regptr[max]));
	}

	for (i = genind[max]; i < dim; i++)
	{
	    Matrix_t *bigmat = NULL;
	    if (side == 'l' && big)
		bigmat = bigform(stdbas[i], Ngen, op_table, d);
	    for (next = 0; next <= max; next++)
	    {
		if (side == 'r')
		{
		    stdbas[dim] = MatDup(stdbas[i]);
		    mat = stdbas[dim];
		    MatMul(mat, gens[next]);
		}
		else /* side == 'l' -- multiple from left */
		{
		    stdbas[dim] = MatDup(stdbas[genind[next]]);
		    if (big)
			MatMul(stdbas[dim], bigmat);
		    else
			MatMul(stdbas[dim], stdbas[i]);
		    mat = stdbas[dim];
		}
		if (independent(stdbas, mat, dim, piv_table, table, regptr[next]))
		{
		    if (big)
		    {
			bpiv[0][0] = piv_table[dim][0];
			bpiv[0][1] = piv_table[dim][1];
			c = 0;
			for (j = 0; j < n; j++)
			{
			    if (basis[j] == NULL)
				continue;
			    a = independent(&(stdbas[dim]), basis[j], 1, bpiv, NULL, NULL);
			    if (c) continue;
			    b = independent(basis, basis[j], j, baspiv, NULL, NULL);
			    if (!a || !b)
			    {
				MatFree(basis[j]);
				basis[j] = NULL;
				c = 1;
			    }
			}
		    }
		    dim++;
    		    MESSAGE(1,("ringgens(): new element3, dim=%d\n",dim));
		}
		else
		    MatFree(stdbas[dim]);

		FfSetNoc(n);
		FfStepPtr(&(regptr[next]));
	    }
	    if (big && side == 'l')
		MatFree(bigmat);
	}
	max++;
    }

    if (!big)
    {
	for (j = 0; j < n; j++)
	    MatFree(basis[j]);
    }

    if (side == 'l')	/* left repr. must be transposed */
    {
	for (i = 0; i < max; i++)
	{
	    mat = MatTransposed(regrep[i]);
	    MatFree(regrep[i]);
	    regrep[i] = mat;
	}
    }
    regrep[max] = NULL;
    gens[max] = NULL;
    return (gens);
}



/* ========================================================================== */






/* ========================================================================== */





static int ParseArgs()

{
    int tmp;

    /* Options.
       -------- */
    standard = AppGetOption(App,"-t");
    hominstd = AppGetOption(App,"-s");
    tmp = AppGetIntOption(App,"-r",0,1,2);
    if (tmp != 0)
    {
	hominstd = 1;
	standard = 1;
	reg = 1;
	side = (char )(tmp == 1 ? 'l' : 'r');
    }
    big = AppGetIntOption(App,"-b",0,0,2);
    if (big != 0)
    {
	if ((TempDirName = AppCreateTempDir(App)) == NULL)
	    return -1;
    }
    hd = AppGetIntOption(App,"-H",0,1,1000000);


    /* Arguments.
       ---------- */
    if (AppGetArguments(App,3,3) < 0)
	return -1;
    MName = App->ArgV[0];
    NName = App->ArgV[1];
    HomName = App->ArgV[2];
    comp = strcmp(MName,NName);
    if (hominstd && comp)
    {
	MTX_ERROR("-b requires <M> = <N>");
	return -1;
    }

    return 0;
}



static int ReadFiles()

{
    /* Read the .cfinfo files.
       ----------------------- */
    if (Lat_ReadInfo(&MInfo,MName) != 0)
	return -1;

    /* Load the generators.
       -------------------- */
    MESSAGE(1,("Reading generators\n"));
    if ((MRep = MrLoad(MInfo.BaseName,MInfo.NGen)) == NULL)
	return -1;
    dim = FfNoc;
    if (comp)
    {
        if ((NRep = MrLoad(NName,MInfo.NGen)) == NULL)
	    return -1;
    }
    else
	NRep = MRep;

    /* Read the head, if -H is used.
       ----------------------------- */
    if (hd > 0)
    {
	Matrix_t *tmp;
	char fn[200];
	sprintf(fn, "%s.rad",MName);
	MESSAGE(1,("Reading the head (%s)\n",fn));
	if ((tmp = MatLoad(fn)) == NULL)
	    return -1;
	if ((rad = MatCutRows(tmp,hd,dim-hd)) == NULL)
	    return -1;
	MatFree(tmp);
	MatEchelonize(rad);
    }

    return 0;
}





static int AllocateWorkspace()

{
    int i;

    FfSetNoc(dim);

    if (   (basis = FfAlloc(dim + 1)) == NULL
        || (space = FfAlloc(dim + 1)) == NULL)
	return -1;
    if ((piv = NALLOC(int,dim + 2)) == NULL)
	return -1;
    if ((op = NALLOC(long,2*dim + 2)) == NULL)
	return -1;
    if (   (stdgen = NALLOC(PTR,MInfo.NGen)) == NULL 
        || (stdtab = NALLOC(long *,MInfo.NGen)) == NULL)
	return -1;
    for (i = 0; i < MInfo.NGen; i++)
    {
	if (   (stdgen[i] = FfAlloc(0)) == NULL
	    || (stdtab[i] = ALLOC(long)) == NULL)
	    return -1;
	stdtab[i][0] = 0;
    }
    if (big == 0)
    {
	if ((tab = ALLOC(long)) == NULL)
	    return -1;
	tab[0] = 0;
    }
    if ((oldstdbas = ALLOC(Matrix_t *)) == NULL)
	return -1;
    if ((oldstdbas[0] = MatAlloc(FfOrder,0,NRep->Gen[0]->Noc)) == NULL)
	return -1;

    return 0;
}


static int Init(int argc, const char **argv)

{
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (ParseArgs() != 0)
	return -1;
    if (ReadFiles() != 0)
	return -1;
    if (AllocateWorkspace() != 0)
    {
	MTX_ERROR("Cannot allocate work space");
	return -1;
    }
    return 0;
}


static void Cleanup()

{
    if (MRep != NULL)
	MrFree(MRep);
    if (NRep != NULL && NRep != MRep)
	MrFree(NRep);
    AppFree(App);
}







/* =========================================================================
   spinpartstdbas() spin up a newdim dimensional part of the 
   standard basis generated by vec beginning at partdim
   ========================================================================= */

Matrix_t *spinpartstdbas(PTR vec, const long *op_table, Matrix_t *gens[], 
    int part_dim, int newdim)

{
    Matrix_t *mat;
    PTR ptr, row;
    int l, newpartdim;

    newpartdim = newdim + part_dim;
    if ((mat = MatAlloc(FfOrder,newdim,gens[0]->Noc)) == NULL)
    {
	MTX_ERROR("Cannot allocate workspace");
	return NULL;
    }
    ptr = mat->Data;
    FfCopyRow(ptr,vec);
    FfStepPtr(&ptr);
    for (l = part_dim + 2; l <= newpartdim; l++)
    {
	row = MatGetPtr(mat,op_table[2*l] - 1 - part_dim);
	FfMapRow(row, gens[op_table[2*l + 1] - 1]->Data, gens[0]->Nor, ptr);
	FfStepPtr(&ptr);
    }
    return mat;
}







/* =========================================================================
   veccont() checks if vec is contained in the subspace generated by mat
   ========================================================================= */

int veccont(Matrix_t *mat, PTR vec, int *pivot_table)

{
    PTR v;
    FEL f;
    int is_contained;

    FfSetNoc(mat->Noc);
    v = FfAlloc(1);
    FfCopyRow(v,vec);
    FfCleanRow(v,mat->Data,mat->Nor,pivot_table);
    is_contained = FfFindPivot(v,&f) < 0;
    SysFree(v);
    return is_contained;
}




static int MakeKernels(int cf, Matrix_t **ker1, Matrix_t **ker2)

{
    char file_name[200];

    /* Load the peak word kernel for M.
       -------------------------------- */
    sprintf(file_name, "%s%s.k", MName, Lat_CfName(&MInfo,cf));
    if ((*ker1 = MatLoad(file_name)) == NULL)
    {
	MTX_ERROR2("Cannot load %s -- did you run 'pwkond %s'?",
	    file_name,MName);
	return 1;
    }

    /* If N is different from M, find the stable peak word kernel in N.
       ---------------------------------------------------------------- */
    if (!comp)
	*ker2 = *ker1;
    else
    {
	Matrix_t *word2;
	WgData_t *wg = WgAlloc(NRep);
	MESSAGE(1,("Calculating the stable peak word kernel in %s\n",NName));
	word2 = WgMakeWord(wg,MInfo.Cf[cf].peakword);
	WgFree(wg);
	MatInsert_(word2, MInfo.Cf[cf].peakpol);

	StablePower_(word2,NULL,ker2);
	MatFree(word2);
    }

    return 0;
}


/* ------------------------------------------------------------------------ */

int main(int argc, const char **argv)

{
    int rc = 0;		/* Program exit code */
    Matrix_t 
	     *homom,	/* the standard basis */
	     *ker1,	/* nullspace in M */
	     *ker2,	/* nullspace in N */
	     **currstdbas,
	     **tempbas,
	     *result,
	     *esys,
	     **gens,
	     **stdbas = NULL,
	     *ech,
	     *echker = NULL;
    PTR basptr, kerptr, echkerptr,
	row, stdgenptr, oldptr, sysptr, resptr, trptr, echptr;
    int hom, sb, seedcount = 0, newdim,
	 newpartdim, homdim = 0, col, size;
    char name[100];
    int i, l;
    int *kerpiv = NULL, *echpiv;
    FEL f;


    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return 1;
    }


    /* Main loop: for each constituent of M.
       ------------------------------------- */
    for (i = 0; i < MInfo.NCf; i++)
    {
	int j;

	MESSAGE(0,("Next constituent: %s%s\n",MName, Lat_CfName(&MInfo,i)));

	if (MakeKernels(i,&ker1,&ker2) != 0)
	    return 1;

	seedcount = 0;
	if (hd)
	{
	    kerpiv = NALLOC(int,ker1->Nor + 1);
	    echker = MatDup(ker1);
	    echkerptr = echker->Data;
	}

	/* Make the next part of the standard basis in M.
	   ---------------------------------------------- */
	for (j = 0; j < ker1->Nor; j++)
	{
	    int k;
	    seedcount++;
	    MESSAGE(1,("Taking kernel vector %d\n",j+1));
	    FfSetNoc(dim);
	    if (hd)
	    {
		PTR ptr;
		FfCleanRow(echkerptr, rad->Data, rad->Nor, rad->PivotTable);
		for (k = 0, ptr = echker->Data; k < j; k++, FfStepPtr(&ptr))
		{
		    if (   kerpiv[k] >= 0 
			&& (f = FfExtract(echkerptr, kerpiv[k])) != FF_ZERO)
		    {
			f = FfDiv(f, FfExtract(ptr, kerpiv[k]));
			FfAddMulRow(echkerptr, ptr, FfNeg(f));
		    }
		}
		if ((kerpiv[j] = FfFindPivot(echkerptr, &f)) < 0)
		{
		    FfStepPtr(&echkerptr);
		    continue;
		}
		FfStepPtr(&echkerptr);
	    }
	    if ((newpartdim = zgensbasis(ker1->Data, seedcount, MInfo.NGen, MRep->Gen, 
		space, piv, basis, partdim, op, stdgen, stdtab)) == partdim) 
	    {
		MESSAGE(1,("No new basis vectors - skipping\n"));
		continue;
	    }
	    MESSAGE(1,("Vector %d (seedcount=%d) spins up to %d\n",j+1,
		seedcount,newpartdim));
	    newdim = newpartdim - partdim;
	    if (!big)
	    {
		if ((tab = NREALLOC(tab,long,++tab[0] + 1)) == NULL)
		    return 1;
		tab[tab[0]] = newdim;
	    }


/* -------------------------------------------------
   extending the standard basis in the second module
   ------------------------------------------------- */
	    MESSAGE(1,("Calculating the standard basis in %s\n",NName));
	    if ((currstdbas = NALLOC(Matrix_t *,ker2->Nor)) == NULL)
		return 1;
	    switch (big)
	    {
		case 0:
		    kerptr = ker2->Data;
		    for (k = 0; k < ker2->Nor; k++)
		    {
			currstdbas[k] = spinpartstdbas(kerptr, op, NRep->Gen, partdim, newdim);
			FfStepPtr(&kerptr);
		    }
		    break;
		case 1:
		    break;
		case 2:
		    kerptr = ker2->Data;
		    for (k = 0; k < ker2->Nor; k++)
		    {
			Matrix_t *mat = spinpartstdbas(kerptr, op, NRep->Gen, partdim, newdim);
			sprintf(name, "%s/curr.%d", TempDirName, k);
			MatSave(mat, name);
			MatFree(mat);
			FfStepPtr(&kerptr);
		    }
	    }


/* ------------------------------------------------------------------------
   builds up the system of equations in order to find the new homomorphisms
   ------------------------------------------------------------------------ */
	    if ((esys = MatAlloc(FfOrder, homdim + ker2->Nor, 
		NRep->Gen[0]->Noc)) == NULL)
		return 1;
	    MESSAGE(1,("Building equation system (%dx%d)\n",esys->Noc,
		esys->Nor));

	    if (esys->Nor == 0)   /* there are no homomorphisms */
	    {
		if (newpartdim == MRep->Gen[0]->Nor)
		{
		    MESSAGE(0,("Warning: There are no homomorphisms from "
			"%s to %s\n",MName,NName));
		   return 0;
		}
		partdim = newpartdim;
		SysFree(currstdbas);
		for (k = 0; k < MInfo.NGen; k++)
		{
		    SysFree(stdgen[k]);
		    if ((stdgen[k] = FfAlloc(0)) == NULL)
			return 1;
		    stdtab[k][0] = 0;
		}
		continue;
	    }

	    if ((ech = MatAlloc(FfOrder, esys->Nor, esys->Nor)) == NULL)
		return 1;
	    echptr = ech->Data;
	    if ((echpiv = NALLOC(int,ech->Nor + 2)) == NULL)
		return 1;
	    ech->Nor = 0;
	    for (k = 0; k < MInfo.NGen; k++)	/* loop for the generators */
	    {
		stdgenptr = stdgen[k];
					/* loop for the operation */

/* the equations for one vector
   ---------------------------- */
		for (l = 1; l <= stdtab[k][0]; l++)
		{
		    int t;
		    Matrix_t *tresys;

		    FfSetNoc(NRep->Gen[0]->Noc);
		    sysptr = esys->Data;

/* the part of the vector in the old submodule
   ------------------------------------------- */
		    for (hom = 0; hom < homdim; hom++) /* the old homomorphisms */
		    {
			Matrix_t *mat = NULL;
			FfMulRow(sysptr, FF_ZERO);
			switch (big)
			{
			    case 0:
				oldptr = oldstdbas[hom]->Data;
				break;
			    case 1:
				mat = bigform(oldstdbas[hom], NRep->Gen, op, partdim);
				oldptr = mat->Data;
				break;
			    case 2:
				sprintf(name, "%s/old.%d", TempDirName, hom);
				mat = MatLoad(name);
				oldptr = mat->Data;
				break;
			}
			for (sb = 0; sb < partdim; sb++)
			{
			    if ((f = FfExtract(stdgenptr, sb)) != FF_ZERO)
				FfAddMulRow(sysptr, oldptr, f);
			    FfStepPtr(&oldptr);
			}
			FfStepPtr(&sysptr);
			if (big)
			    MatFree(mat);
		    }

/* the part of the vector over the old submodule
   --------------------------------------------- */
		    kerptr = ker2->Data;
		    for (hom = 0; hom < ker2->Nor; hom++)
		    {
			switch (big)
			{
			    case 0:
				break;
			    case 1:
				currstdbas[hom] = spinpartstdbas(kerptr, op, NRep->Gen, partdim, newdim);
				break;
			    case 2:
				sprintf(name, "%s/curr.%d", TempDirName, hom);
				currstdbas[hom] = MatLoad(name);
				break;
			}
			basptr = MatGetPtr(currstdbas[hom],stdtab[k][l] - 1);
			FfMapRow(basptr, NRep->Gen[k]->Data, NRep->Gen[0]->Nor, sysptr);
			FfMulRow(sysptr, FfNeg(FF_ONE));
			basptr = currstdbas[hom]->Data;
			for (sb = partdim; sb < newpartdim; sb++)
			{
			    if ((f = FfExtract(stdgenptr, sb)) != FF_ZERO)
				FfAddMulRow(sysptr, basptr, f);
			    FfStepPtr(&basptr);
			}
			FfStepPtr(&sysptr);
			switch (big)
			{
			    case 1:
				FfStepPtr(&kerptr);
				MatFree(currstdbas[hom]);
				break;
			    case 2:
				MatFree(currstdbas[hom]);
				break;
			}
		    }
		    FfSetNoc(dim);
		    FfStepPtr(&stdgenptr);

/* ------------------------------------
   eliminating the superfluous equations
   ------------------------------------ */
		    tresys =  MatTransposed(esys);
		    trptr = tresys->Data;
		    FfSetNoc(tresys->Noc);
		    for (t = 0; t < tresys->Nor; t++)
		    {
			FfCleanRow(trptr, ech->Data , ech->Nor, echpiv);

			if ((echpiv[ech->Nor] = FfFindPivot(trptr, &f)) >= 0)
			{
			    FfCopyRow(echptr, trptr);
			    if ((++ech->Nor) > ech->Noc) 
			    {
				MTX_ERROR("The matrix has rank greater than number of rows");
				return 1;
			    }
			    FfStepPtr(&echptr);
			}
			FfStepPtr(&trptr);
		    }
		    MatFree(tresys);
		}	/* end of loop with l */
	    }		/* end of loop with k */


/* ----------------------------------------
   solving the remained system of equations
   ---------------------------------------- */
	    if (big == 2)
	    {
		for (k = 0; k < homdim; k++)
		{
		    sprintf(name, "%s/old.%d", TempDirName, k);
		    if (remove(name))
			printf("Error by removing file %s\n", name);
		}
	    }
	    if ((ech->Data = (PTR) SysRealloc(ech->Data,FfCurrentRowSize*ech->Nor)) == NULL)
		return 1;
	    MESSAGE(1,("Solving equation system (%dx%d)\n",ech->Nor,ech->Noc));
	    if (ech->Nor > 0)
		result = MatNullSpace__(MatTransposed(ech));
	    else
		result = MatId(FfOrder, ech->Noc);
	    MatFree(ech);
	    MatFree(esys);
	    tempbas = oldstdbas;
	    if (   (oldstdbas = NALLOC(Matrix_t *,result->Nor)) == NULL
		|| (oldstdbas[0] = MatAlloc(FfOrder,0, NRep->Gen[0]->Noc)) == NULL)
		return 1;

/* --------------------------------------
   extending the extendable homomorphisms
   -------------------------------------- */

	    size = (big == 0 ? newpartdim : tempbas[0]->Nor + 1);
	    resptr = result->Data;
	    for (k = 0; k < result->Nor; k++)
	    {
		FfSetNoc(NRep->Gen[0]->Noc);
		if ((oldstdbas[k] = MatAlloc(FfOrder, size, NRep->Gen[0]->Noc)) == NULL)
		    return 1;
		for (l = 0; l < homdim; l++)
		{
		    int t;
		    int oldnoc = FfNoc;
		    oldptr = oldstdbas[k]->Data;
		    row = tempbas[l]->Data;
		    FfSetNoc(result->Noc);
		    f = FfExtract(resptr,l);
		    FfSetNoc(oldnoc);
		    for (t = 0; t < tempbas[0]->Nor; t++)
		    {
			FfAddMulRow(oldptr, row, f);
			FfStepPtr(&oldptr);
			FfStepPtr(&row);
		    }
		}
		basptr = MatGetPtr(oldstdbas[k],tempbas[0]->Nor);
		if (big)
		    row = ker2->Data;
		for (l = 0, col = homdim; l < ker2->Nor; l++, col++)
		{
		    int t;
		    int oldnoc = FfNoc;
		    oldptr = basptr;
		    FfSetNoc(result->Noc);
		    f = FfExtract(resptr, col);
		    FfSetNoc(oldnoc);
		    switch (big)
		    {
			case 0:
			    row = currstdbas[l]->Data;
			    for (t = 0; t < newdim; t++)
			    {
				FfAddMulRow(oldptr, row, f);
				FfStepPtr(&oldptr);
				FfStepPtr(&row);
			    }
			    break;
			case 1:
			case 2:
			    FfAddMulRow(oldptr, row, f);
			    FfStepPtr(&row);
			    break;
		    }
		}

		FfSetNoc(result->Noc);
		FfStepPtr(&resptr);
		if (big == 2 && newpartdim < dim)
		{
		    Matrix_t *mat = bigform(oldstdbas[k], NRep->Gen, op, newpartdim);
		    sprintf(name, "%s/old.%d", TempDirName, k);
		    MatSave(mat, name);
		    MatFree(mat);
		}
	    }

/* gives back the superfluous space to the memory 
   ---------------------------------------------- */
	    for (k = 0; k < homdim; k++)
		MatFree(tempbas[k]);
	    SysFree(tempbas);
	    switch (big)
	    {
		case 0:
		    for (k = 0; k < ker2->Nor; k++)
			MatFree(currstdbas[k]);
		    SysFree(currstdbas);
		    break;
		case 1:
		case 2:
		    SysFree(currstdbas);
		    break;
	    }
	    homdim = result->Nor;
	    MESSAGE(0,("%d homomorphisms found\n",homdim));
	    MatFree(result);
	    partdim = newpartdim;
	    for (k = 0; k < MInfo.NGen; k++)
	    {
		SysFree(stdgen[k]);
		if ((stdgen[k] = FfAlloc(0)) == NULL)
		    return 1;
		stdtab[k][0] = 0;
	    }


	    if (newpartdim == dim)
	    {
/* ---------------
   makes the output
   --------------- */

		homom = MatAlloc(FfOrder,dim,dim);
		SysFree(homom->Data);
		homom->Data = basis;
		sprintf(name,"%s.std",MName);
		MESSAGE(1,("Writing standard basis to %s\n",name));
		MatSave(homom,name);

		if (standard || hominstd)
		{
		    Matrix_t *homomi = MatInverse(homom);
		    if (standard)
		    {
			MESSAGE(1,("Transforming %s into standard basis\n",
			    MName));
			for (k = 0; k < MInfo.NGen; k++)
			{
			    Matrix_t *mat = MatDup(homom);
			    MatMul(mat, MRep->Gen[k]);
			    MatMul(mat, homomi);
			    sprintf(name, "%s.std.%d", MName, k + 1);
			    MatSave(mat, name);
			    if (reg)
			    {
				MatFree(NRep->Gen[k]);
				NRep->Gen[k] = mat;
			    }
			}
		    }
		    if (hominstd)
		    {
			MESSAGE(1,("Transforming homomorphisms into standard "
			    "basis\n"));
			for (k = 0; k < homdim; k++)
			    MatMul(oldstdbas[k], homomi);
		    }
		}


		if (reg)
		{
		    Matrix_t **regrep;
		    Lat_Info end_info;

		    MESSAGE(1,("Calculating regular representation\n"));
		    if ((regrep = NALLOC(Matrix_t *,homdim)) == NULL)
			return 1;
		    if ((stdbas = NALLOC(Matrix_t *,homdim + 1)) == NULL)
			return 1;
		    gens = ringgens(oldstdbas,homdim,tab,regrep,side,big,stdbas,op,NRep->Gen);

/* gens werden voellig gemacht. vielleicht kann es noch geaendert werden. */


		    for (k = 0; gens[k] != NULL; k++)
		    {
			sprintf(name, "%s.gens.%d", HomName, k + 1);
			MatSave(gens[k], name);
			sprintf(name, "%s.%crr.%d", HomName, side, k + 1);
			MatSave(regrep[k], name);
		    }

		    /* Cretae the <endo>.lrr.cfinfo file */
		    memset(&end_info,0,sizeof(end_info));
		    end_info.NGen = k;
		    sprintf(end_info.BaseName,"%s.%crr",HomName,side);
		    Lat_WriteInfo(&end_info);
		}
		MESSAGE(1,("Writing homomorphisms\n"));
		for (k = 0; k < homdim; k++)
		{
		    Matrix_t *mat = (reg ? stdbas[k] : oldstdbas[k]);
		    sprintf(name, "%s.%d", HomName, k + 1);
		    switch (big)
		    {
			case 0:
			    MatSave(mat, name);
			    break;
			case 1:
			case 2:
			    {
				Matrix_t *m = bigform(mat, NRep->Gen, op, NRep->Gen[0]->Noc);
				MatSave(m, name);
				MatFree(m);
			    }
			    break;
		    }
		}
		if (big == 2)
		{
		    for (k = 0; k < ker2->Nor; k++)
		    {
			sprintf(name, "%s/curr.%d", TempDirName, k);
			SysRemoveFile(name);
		    }
		}
		Cleanup();
		return 0;
	    }
	    if (big == 2)
	    {
		for (k = 0; k < ker2->Nor; k++)
		{
		    sprintf(name, "%s/curr.%d", TempDirName, k);
		    SysRemoveFile(name);
		}
	    }
	}		/* end of loop with j */
	MatFree(ker2);
	if (comp) 
	    MatFree(ker1);
    }


    Cleanup();
    return rc;
}



/**
!section apps.lattice
!program mkhom "Homomorphisms"
!synopsis 
    mkhom [-ts] [-r <Side>] [-b <Mode>] [-H <Dim>] <M> <N> <Hom>
!option -t
    Calculate standard generators for {\it M}.
!option -s
    When <M> equals <N>, give endomorphisms in standard basis.
!option -r
    When <M> equals <N>, find a generating set of $\End(M)$, and
    calculate the left ({\it Side}=1) or right ({\it Side}=2) regular 
    representation.
!option -b
    Save memory. <Mode> can be $0$, $1$, or $2$ (see below).
!option -H
    If the radical is given, <Dim> is the dimension of the head.
 ** @param M
    First representation.
 ** @param N
    Second representation.
 ** @param Hom
    Output file name.
 ** @see chop pwkond
!description
    This program calculates a basis for the vector space of homomorphisms
    between two $kG$-modules, $\Hom_{kG}(M,N)$. In the case $M=N$ the program
    optionally finds a generating set for the algebra of endomoprhisms,
    $\End_{kG}(M)$, and calculates the corresponding left or right regular 
    representation. 
    
    If used without any options, MKHOM writes the standard basis of $M$
    to `<M>.std', and a $k$-basis of the homomorphism space to 
    `<Hom>.1', `<Hom>.2', ... The latter are given with 
    respect to the standard basis of $M$ and the original basis of $N$.
    To get the homomorphisms with respect to the original bases of $M$
    and $N$, multiply the matrices from the left with the inverse of 
    `<Hom>'.

    MKHOM uses peak words of the first module. Thus, before using the program, 
    CHOP and PWKOND must have been run on the first module.
    
!implementation
    The algorithm used by this program was developed by Magdolna Sz"oke
<<<<<<< HEAD
    \cite{Sz98}.
=======
    @ref Sz98.
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
 **/



