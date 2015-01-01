////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Calculate radical series or homomorphism of the PIMs 
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>



/* --------------------------------------------------------------------------
   Global data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO

static MtxApplicationInfo_t AppInfo = { 
"rad", "Radical series", 
"SYNTAX\n"
"    rad " MTX_COMMON_OPTIONS_SYNTAX " [-l <Length>] [-H <Num>] <Name>\n"
"\n"
"FILES\n"
"    <Name>.1 ... <Name>.nbgen	i  generators of a representation\n"
"    <Name>.cfinfo		i  info-file after running PWKOND\n"
"    <Name>.rad			o  matrix for basischange\n"
"				   or the vectors for hom\n"
"    <Name><S>.h<Num>		o  ???\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -l <Length> ............. calculate the first <length> heads\n"
"    -H <Num> ................ calculate generators for the <num>th head\n"
"\n"
"DESCRIPTION\n"
"\n"
"    This program calculates the radical series of an arbitrary\n"
"    module given by <gens> or the homomorphisms from the\n"
"    projective modules corresponding to the composition factors\n"
"    of the given module to the module.\n" 
"    In case of option 'h', the vectors generating submodules with\n"
"    head isomorphic to <cfs>, are stored in gens<Name>.h<num>\n"
};

static MtxApplication_t *App = NULL;
static const char *Name = NULL;
static Lat_Info Info;
static int Head = 0;		/* Make head */
static int Len = 0;		/* Max. length of radical series */
static MatRep_t *Rep;		/* The representation (on M) */
static MatRep_t *CfRep[LAT_MAXCF];  /* Constituents in standard basis */
static Matrix_t *sed[LAT_MAXCF];   /* Kernel of the peak words */
static IntMatrix_t *OpTable[LAT_MAXCF];    /* Operations, written to <Op> */



/*-----------------------------------------------------------------------*/




Matrix_t *intersect(Matrix_t *mat1,Matrix_t *mat2)

{
    PTR wrk1, wrk2;
    int nor1 = mat1->Nor, nor2 = mat2->Nor;
    int *piv = NALLOC(int,nor1 + nor2);
    Matrix_t *result;

    FfSetNoc(mat1->Noc);
    wrk1 = FfAlloc(nor1 + nor2);
    wrk2 = FfAlloc(nor1 + nor2);
    memcpy(wrk1,mat1->Data,nor1 * FfCurrentRowSize);
    memcpy(FfGetPtr(wrk1,nor1),mat2->Data,nor2 * FfCurrentRowSize);

    FfSumAndIntersection(wrk1,&nor1,&nor2,wrk2,piv);

    result = MatAlloc(FfOrder,nor2,FfNoc);
    memcpy(result->Data,FfGetPtr(wrk2,nor1),nor2 * FfCurrentRowSize);

    SysFree(wrk1);
    SysFree(wrk2);
    SysFree(piv);
    return result;
}




/* making the dual of a representation given by gens 
   -------------------------------------------------- */

static void Dualize(MatRep_t *rep)

{
    int i;
    for (i=0; i < rep->NGen; i++)
    {
	Matrix_t *mat = MatTransposed(rep->Gen[i]);
	MatFree(rep->Gen[i]);
	rep->Gen[i] = mat;
    }
}




static int ReadFiles()

{
    int i;

    /* Read input files.
       ----------------- */
    Name = App->ArgV[0];
    if (Lat_ReadInfo(&Info,Name) != 0)
	return -1;
    if (!Head)
	Info.NHeads = 0;

    /* Read the generators for the module.
       ----------------------------------- */
    Rep = MrLoad(Name,Info.NGen);
    if (Rep == NULL)
	return -1;

    
    for (i = 0; i < Info.NCf; ++i)
    {
	char fn[200];
	/* Read generators
	   --------------- */
	sprintf(fn,"%s%s.std", Info.BaseName,Lat_CfName(&Info,i));
	CfRep[i] = MrLoad(fn,Info.NGen);
	if (CfRep[i] == NULL)
	    return -1;
    }

    return 0;
}



static int Init(int argc, const char **argv)

{
    int headno, maxlen;

    /* Process command line options.
       ----------------------------- */
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    headno = AppGetIntOption(App,"-H --head",-1,1,1000);
    maxlen = AppGetIntOption(App,"-l --max-length",-1,1,1000);
    if (headno != -1 && maxlen != -1)
    {
	MTX_ERROR("'-l' and '-H' cannot be used together");
	return -1;
    }
    if (headno != -1)
    {
	Head = 1;
	Len = headno;
    }
    else
	Len = maxlen;

    /* Process command line arguments.
       ------------------------------- */
    if (AppGetArguments(App,1,1) != 1)
	return -1;
    if (ReadFiles() != 0)
    {
	MTX_ERROR("Error reading input files");
	return 1;
    }
    return 0;
}


static int DualizeConstituents()

{
    int j;

    for (j = 0; j < Info.NCf; ++j)
    {
	Matrix_t *word, *word_tr, *sb, *seed;

	WgData_t *wg = WgAlloc(CfRep[j]);
	if (Info.Cf[j].peakword <= 0)
	{
	    MTX_ERROR("No peak word found. Run PWKOND first!");
	    return -1;
	}
	word = WgMakeWord(wg,Info.Cf[j].peakword);
	word_tr = MatTransposed(word);
	MatFree(word);
	seed = MatNullSpace__(MatInsert_(word_tr,Info.Cf[j].peakpol));
	WgFree(wg);
	Dualize(CfRep[j]);
	OpTable[j] = NULL;
	sb = SpinUp(seed,CfRep[j],SF_FIRST|SF_CYCLIC|SF_STD,OpTable + j,NULL);
	if (sb == NULL)
	{
	    MTX_ERROR("Cannot make standard basis");
	    return -1;
	}
	MrChangeBasis(CfRep[j],sb);
	MatFree(sb);
    }

    return 0;
}

 
int main(int argc, const char **argv)

{
    int i, j, socdim=0, soclen = 0, dim;
    int flag;
    int *cfvec;
    char name[LAT_MAXBASENAME];
    Matrix_t *stgen, *partbas, 
             *mat, *seed, 
	     *emb = NULL, *basis = NULL,
	     *soc2 = NULL, *rad2, *echbas;
    WgData_t *rep;

    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return 1;
    }
    if (DualizeConstituents() != 0)
    {
	MTX_ERROR("Error while dualizing constituents");
	return 1;
    }


    while(1) 
    {
	Matrix_t *bas, *basi;
	cfvec = NALLOC(int,Info.NCf);
	bas = MatAlloc(FfOrder, Rep->Gen[0]->Nor, Rep->Gen[0]->Noc);
	if (cfvec == NULL || bas == NULL)
	{
	    MTX_ERROR("Cannot allocate work space");
	    return 1;
	}

/* -------------------------------------------------------------
   determining the nullspace of the peakwords in the dual module
   ------------------------------------------------------------- */
	rep = WgAlloc(Rep);
	for (j = 0; j < Info.NCf; j++)
	{
	    Matrix_t *word, *w;
	    word = WgMakeWord(rep,Info.Cf[j].peakword);
	    w = MatTransposed(word);
	    MatFree(word);
            sed[j] = MatNullSpace__(MatInsert_(w, Info.Cf[j].peakpol)); 
	}
	WgFree(rep);


	Dualize(Rep);

	for (j=0; j < Info.NCf; j++)
	{

/* ------------------------------------------------------------------
   computes the submodules isomorphic to the given composition factor
   ------------------------------------------------------------------ */

            if (sed[j]->Nor != 0) 
	    {
                partbas = HomogeneousPart(Rep,CfRep[j],sed[j],OpTable[j],
		    Info.Cf[j].spl);
		MatFree(sed[j]);
	    }
            else 
		partbas = sed[j];

	    cfvec[j] = partbas->Nor / Info.Cf[j].dim;
	    MatCopyRegion(bas,socdim,0,partbas,0,0,-1,-1);
            socdim += partbas->Nor;
	    MESSAGE(2,("  headdim of the first %d cfs is %d\n", j+1, socdim));
            MatFree(partbas);
        }

/* makes the output
   ---------------- */
	++soclen;
	MESSAGE(0,("Head %d: %d =",soclen,socdim));
	flag = 0;
	for (j = 0; j < Info.NCf; j++)
	{
	    if (cfvec[j] <= 0) 
		continue;
	    if (flag++ > 0)
		MESSAGE(0,(" +"));
	    if (cfvec[j] == 1)
		MESSAGE(0,(" %s",Lat_CfName(&Info,j)));
	    else
		MESSAGE(0,(" %d*%s",cfvec[j],Lat_CfName(&Info,j)));
	}
	MESSAGE(0,("\n"));
	if (!Head)
	    Lat_AddHead(&Info,cfvec);
	SysFree(cfvec);


/* ------------------------------------------------------------
   makes the socle in the factormodule and leaves in case of -H
   ------------------------------------------------------------ */

	if (Head && Len == soclen)
	{
	    soc2 = MatCut(bas,0,0,socdim,bas->Noc);
	    MatFree(bas);
	    break;
	}


/* -----------------------------------
   exiting if the module is semisimple
   ----------------------------------- */
	if (socdim == Rep->Gen[0]->Nor) 
        {
	    stgen = MatInverse(bas);
	    MatFree(bas);
	    bas = MatTransposed(stgen);
	    MatFree(stgen);
	    sprintf(name,"%s.rad",Name);
	    if (basis == NULL)
		MatSave(bas,name);
	    else
	    {
		Matrix_t *mat = MatCutRows(basis,basis->Nor - socdim,socdim);
		MatMul(bas,mat);
		MatFree(mat);
		MatCopyRegion(basis,basis->Nor - socdim,0,bas,0,0,-1,-1);
		MatSave(basis, name);
	    }
            break;
        }

/* -------------------------------------------------------------
   extends the basis of the socle to a basis of the whole module
   ------------------------------------------------------------- */

	MatEchelonize(bas);
	echbas = MatAlloc(bas->Field,bas->Noc,bas->Noc);
	MatCopyRegion(echbas,0,0,bas,0,0,-1,-1);
	for (i = bas->Nor; i < bas->Noc; ++i)
	    FfInsert(MatGetPtr(echbas,i),bas->PivotTable[i],FF_ONE);
	MatFree(bas);
	bas = echbas;


        basi = MatInverse(bas);
	dim = bas->Nor - socdim;
	stgen = MatTransposed(basi);

/* multiplying the last two basis-changes
   -------------------------------------- */

	/*if(1)	*************************************************************/
	{
	    if(basis == NULL)
		basis = stgen;
	    else
	    {
		Matrix_t *mat;
		mat = MatCutRows(basis,basis->Nor - stgen->Nor,stgen->Nor);
		MatMul(stgen, mat);
		MatCopyRegion(basis,basis->Nor - stgen->Nor,0,stgen,0,0,-1,-1);
		MatFree(mat);
		MatFree(stgen);
	    }
	}

	if (Len == soclen && !Head)
	{
	    sprintf(name, "%s.rad", Name);
	    MatSave(basis, name);
	    break;
	}

/* calculates the embedding in the Len-1st radical, in case of -H
   -------------------------------------------------------------- */
	if (Len == soclen + 1 && Head)
	{
	    emb = MatCutRows(basis,basis->Nor - dim,dim);
	    if (emb == NULL)
	    {
		MTX_ERROR("Cannot allocate matrxi");
		return 1;
	    }
	}

/* --------------------------------------
   basis transformation and factorization
   -------------------------------------- */

	MESSAGE(1,("Reducing to dimension %d\n",Rep->Gen[0]->Noc-socdim));
        for (i = 0; i < Info.NGen; ++i)
        {
            stgen = MatDup(bas);
            MatMul(stgen,Rep->Gen[i]);
            MatMul(stgen,basi);		    /* the transformation */
            MatFree(Rep->Gen[i]);

	    partbas = MatCutRows(stgen,socdim,dim);
	    MatFree(stgen);
	    stgen = MatTransposed(partbas);	/* 'dualizing' */
	    MatFree(partbas);
	    Rep->Gen[i] = MatCutRows(stgen,socdim,dim);
	    MatFree(stgen);
        }

	MatFree(bas);
        MatFree(basi);
        FfSetNoc(Rep->Gen[0]->Noc);
	socdim = 0;
    }


    Lat_WriteInfo(&Info);



    if (socdim == Rep->Gen[0]->Nor && !Head)
	return 0;	
    if (socdim < Rep->Gen[0]->Nor && !Head) 
    {
	MESSAGE(0,("Radical length is greater than %d\n",soclen));
	if (!Head)
	    return 0;	
     }

    

/*----------------------------------------------------------------------------*/

    if (soc2 == NULL)
    {
	MESSAGE(0,("Radical length is smaller than %d, there are no vectors "
		"in the %dth Head\n",Len,Len));
	return 0;
    }

/* -----------------------------
   computing the <Len>th radical
   ----------------------------- */

    rad2 = MatNullSpace__(MatTransposed(soc2));	/* rad^Len given in rad^(Len - 1) */

    Dualize(Rep);	/* the action on rad^(Len - 1) */

/* ------------------------------------------------------------------
   calculates the vecors generating the irreducibles lying in rad^Len
   ------------------------------------------------------------------ */
    rep = WgAlloc(Rep);
    for(j = 0; j < Info.NCf; j++)
    {
	Matrix_t *word = WgMakeWord(rep,Info.Cf[j].peakword);
	Matrix_t *seed2;

	/* makes the iterated nullspace of the peakword */
	MatInsert_(word, Info.Cf[j].peakpol);
	if ((seed = MatAlloc(FfOrder, 0, 0)) == NULL)
	    return 1;
	for (seed2 = MatNullSpace(word); seed->Nor < seed2->Nor; 
	    seed2 = MatNullSpace(MatMul(word, word)))
	{
	    MatFree(seed);
	    seed = seed2;
	}
	MatFree(seed2);
	MatFree(word);

	if (seed->Nor > 0)
	{
	    Matrix_t *sec = intersect(seed, rad2); /* the nullspace in rad^2 */
	    MatPivotize(sec);
	    MatClean(seed,sec);
	    MatFree(sec);
	    if (emb != NULL)
		MatMul(seed, emb); /* embedding into the original module */
	    mat = seed;
	}
	else
	{
	    int noc = (emb == NULL) ? Rep->Gen[0]->Noc : emb->Noc;
	    MatFree(seed);
	    if ((mat = MatAlloc(FfOrder, 0, noc)) == NULL)
	    {
		MTX_ERROR("Cannot allocate matrix");
		return 1;
	    }
	}

	sprintf(name,"%s%s.h%d", Name, Lat_CfName(&Info,j), Len);
	MatSave(mat,name);
/*	MatPrint(name, mat);*/
	MatFree(mat);
    }
    WgFree(rep);
    return(0);
}



/**
@page prog_rad rad - Radical Series

@section rad_syntax Command Line
<pre>
rad @em Options [-l @em MaxLength] [-H @em Num] @em Module
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -l @em MaxLength
  Maximum number of layers to compute.
@par -H @em Num
  Calculate generators for the @em Num-th head.
@par @em Module
  Module name.

@section rad_inp Input Files
@par @em Name.cfinfo
  Constituent information.
@par @em NameCf.std.1, @em NameCf.std.2, ...
  Generators on the irreducible constituents.

@section rad_out Output Files
@par @em Name.cfinfo
  Radical information.
@par @em NameCf.hX
  Generators for the X-th head.

@section rad_desc Description
This program calculates the radical series of an arbitrary module @em Name,
or the homomorphisms from the projective modules
corresponding to the composition factors of the given module to the module.

@section rad_impl Implementation Details
The program uses an algorithm by Magdolna Szöke, see @ref Sz98 "[Sz98]".

**/
