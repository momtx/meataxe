////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Calculate radical series or homomorphism of the PIMs 
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>



/* --------------------------------------------------------------------------
   Global data
   -------------------------------------------------------------------------- */


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
static LatInfo_t* Info;
static int Head = 0;                    // Make head
static int Len = 0;                     // Max. length of radical series
static MatRep_t *Rep;                   // The representation (on M)
static MatRep_t *CfRep[LAT_MAXCF];      // Constituents in standard basis
static Matrix_t *sed[LAT_MAXCF];        // Kernel of the peak words
static IntMatrix_t *OpTable[LAT_MAXCF]; // Operations, written to <Op> 

////////////////////////////////////////////////////////////////////////////////////////////////////

Matrix_t *intersect(Matrix_t *mat1, Matrix_t *mat2)
{
    PTR wrk1, wrk2;
    uint32_t nor1 = mat1->nor, nor2 = mat2->nor;
    uint32_t *piv = NALLOC(uint32_t,nor1 + nor2);
    Matrix_t *result;

    wrk1 = ffAlloc(nor1 + nor2, mat1->noc);
    wrk2 = ffAlloc(nor1 + nor2, mat1->noc);
    memcpy(wrk1,mat1->data, ffSize(nor1, mat1->noc));
    memcpy(ffGetPtr(wrk1,nor1,mat1->noc),mat2->data,ffSize(nor2, mat1->noc));

    ffSumAndIntersection(mat1->noc, wrk1,&nor1,&nor2,wrk2,piv);

    result = matAlloc(ffOrder,nor2,mat1->noc);
    memcpy(result->data,ffGetPtr(wrk2,nor1,mat1->noc),ffSize(nor2, mat1->noc));

    sysFree(wrk1);
    sysFree(wrk2);
    sysFree(piv);
    return result;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

// making the dual of a representation given by gens 

static void Dualize(MatRep_t *rep)
{
    int i;
    for (i=0; i < rep->NGen; i++)
    {
	Matrix_t *mat = matTransposed(rep->Gen[i]);
	matFree(rep->Gen[i]);
	rep->Gen[i] = mat;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void ReadFiles()
{
    Name = App->argV[0];
    Info = latLoad(Name);
    if (!Head)
	Info->NHeads = 0;

    // Read the generators for the module and the constituents.
    Rep = mrLoad(Name,Info->NGen);
    for (int i = 0; i < Info->nCf; ++i)
    {
	char fn[200];
	sprintf(fn,"%s%s.std", Info->BaseName,latCfName(Info,i));
	CfRep[i] = mrLoad(fn,Info->NGen);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
    int headno, maxlen;

    /* Process command line options.
       ----------------------------- */
    App = appAlloc(&AppInfo,argc,argv);
    headno = appGetIntOption(App,"-H --head",-1,1,1000);
    maxlen = appGetIntOption(App,"-l --max-length",-1,1,1000);
    if (headno != -1 && maxlen != -1)
	mtxAbort(MTX_HERE,"'-l' and '-H' cannot be used together");
    if (headno != -1)
    {
	Head = 1;
	Len = headno;
    }
    else
	Len = maxlen;

    /* Process command line arguments.
       ------------------------------- */
    appGetArguments(App,1,1);
    ReadFiles();
}


static int DualizeConstituents()

{
    int j;

    for (j = 0; j < Info->nCf; ++j)
    {
	Matrix_t *word, *word_tr, *sb, *seed;

	WgData_t *wg = wgAlloc(CfRep[j]);
	if (Info->Cf[j].peakWord == 0)
	    mtxAbort(MTX_HERE,"No peak word found. Run PWKOND first!");
	word = wgMakeWord(wg,Info->Cf[j].peakWord);
	word_tr = matTransposed(word);
	matFree(word);
	seed = matNullSpace__(matInsert_(word_tr,Info->Cf[j].peakPol));
	wgFree(wg);
	Dualize(CfRep[j]);
	OpTable[j] = NULL;
	sb = SpinUp(seed,CfRep[j],SF_FIRST|SF_CYCLIC|SF_STD,OpTable + j,NULL);
	if (sb == NULL)
	    mtxAbort(MTX_HERE,"Cannot make standard basis");
	mrChangeBasis(CfRep[j],sb);
	matFree(sb);
    }

    return 0;
}

 
int main(int argc, char **argv)

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

    init(argc,argv);
    if (DualizeConstituents() != 0)
	mtxAbort(MTX_HERE,"Error while dualizing constituents");

    while(1) 
    {
	Matrix_t *bas, *basi;
	cfvec = NALLOC(int,Info->nCf);
	bas = matAlloc(ffOrder, Rep->Gen[0]->nor, Rep->Gen[0]->noc);

        // determining the nullspace of the peakwords in the dual module
	rep = wgAlloc(Rep);
	for (j = 0; j < Info->nCf; j++)
	{
	    Matrix_t *word, *w;
	    word = wgMakeWord(rep,Info->Cf[j].peakWord);
	    w = matTransposed(word);
	    matFree(word);
            sed[j] = matNullSpace__(matInsert_(w, Info->Cf[j].peakPol)); 
	}
	wgFree(rep);


	Dualize(Rep);

	for (j=0; j < Info->nCf; j++)
	{

/* ------------------------------------------------------------------
   computes the submodules isomorphic to the given composition factor
   ------------------------------------------------------------------ */

            if (sed[j]->nor != 0) 
	    {
                partbas = HomogeneousPart(Rep,CfRep[j],sed[j],OpTable[j],
		    Info->Cf[j].spl);
		matFree(sed[j]);
	    }
            else 
		partbas = sed[j];

	    cfvec[j] = partbas->nor / Info->Cf[j].dim;
	    matCopyRegion(bas,socdim,0,partbas,0,0,partbas->nor, partbas->noc);
            socdim += partbas->nor;
	    MTX_LOG2("  headdim of the first %d cfs is %d\n", j+1, socdim);
            matFree(partbas);
        }

/* makes the output
   ---------------- */
	++soclen;
	MTX_LOGI("Head %d: %d =",soclen,socdim);
	flag = 0;
	for (j = 0; j < Info->nCf; j++)
	{
	    if (cfvec[j] <= 0) 
		continue;
	    if (flag++ > 0)
		MTX_LOGI(" +");
	    if (cfvec[j] == 1)
		MTX_LOGI(" %s",latCfName(Info,j));
	    else
		MTX_LOGI(" %d*%s",cfvec[j],latCfName(Info,j));
	}
	MTX_LOGI("\n");
	if (!Head)
	    latAddHead(Info,cfvec);
	sysFree(cfvec);


/* ------------------------------------------------------------
   makes the socle in the factormodule and leaves in case of -H
   ------------------------------------------------------------ */

	if (Head && Len == soclen)
	{
	    soc2 = matDupRegion(bas,0,0,socdim,bas->noc);
	    matFree(bas);
	    break;
	}


/* -----------------------------------
   exiting if the module is semisimple
   ----------------------------------- */
	if (socdim == Rep->Gen[0]->nor) 
        {
	    stgen = matInverse(bas);
	    matFree(bas);
	    bas = matTransposed(stgen);
	    matFree(stgen);
	    sprintf(name,"%s.rad",Name);
	    if (basis == NULL)
		matSave(bas,name);
	    else
	    {
		Matrix_t *mat = matDupRows(basis,basis->nor - socdim,socdim);
		matMul(bas,mat);
		matFree(mat);
		matCopyRegion(basis,basis->nor - socdim,0,bas,0,0,bas->nor, bas->noc);
		matSave(basis, name);
	    }
            break;
        }

/* -------------------------------------------------------------
   extends the basis of the socle to a basis of the whole module
   ------------------------------------------------------------- */

	matEchelonize(bas);
	echbas = matAlloc(bas->field,bas->noc,bas->noc);
	matCopyRegion(echbas,0,0,bas,0,0,bas->nor, bas->noc);
	for (i = bas->nor; i < bas->noc; ++i)
	    ffInsert(matGetPtr(echbas,i),bas->pivotTable[i],FF_ONE);
	matFree(bas);
	bas = echbas;


        basi = matInverse(bas);
	dim = bas->nor - socdim;
	stgen = matTransposed(basi);

/* multiplying the last two basis-changes
   -------------------------------------- */

	/*if(1)	*************************************************************/
	{
	    if(basis == NULL)
		basis = stgen;
	    else
	    {
		Matrix_t *mat;
		mat = matDupRows(basis,basis->nor - stgen->nor,stgen->nor);
		matMul(stgen, mat);
		matCopyRegion(basis,basis->nor - stgen->nor,0,stgen,0,0,stgen->nor, stgen->noc);
		matFree(mat);
		matFree(stgen);
	    }
	}

	if (Len == soclen && !Head)
	{
	    sprintf(name, "%s.rad", Name);
	    matSave(basis, name);
	    break;
	}

/* calculates the embedding in the Len-1st radical, in case of -H
   -------------------------------------------------------------- */
	if (Len == soclen + 1 && Head)
	{
	    emb = matDupRows(basis,basis->nor - dim,dim);
	    if (emb == NULL)
	    {
		mtxAbort(MTX_HERE,"Cannot allocate matrxi");
		return 1;
	    }
	}

        // basis transformation and factorization
	MTX_LOGD("Reducing to dimension %lu",(unsigned long)(Rep->Gen[0]->noc-socdim));
        for (i = 0; i < Info->NGen; ++i)
        {
            stgen = matDup(bas);
            matMul(stgen,Rep->Gen[i]);
            matMul(stgen,basi);		    // the transformation
            matFree(Rep->Gen[i]);

	    partbas = matDupRows(stgen,socdim,dim);
	    matFree(stgen);
	    stgen = matTransposed(partbas);	// 'dualizing'
	    matFree(partbas);
	    Rep->Gen[i] = matDupRows(stgen,socdim,dim);
	    matFree(stgen);
        }

	matFree(bas);
        matFree(basi);
	socdim = 0;
    }


    latSave(Info);

    if (socdim == Rep->Gen[0]->nor && !Head)
	return 0;	
    if (socdim < Rep->Gen[0]->nor && !Head) 
    {
	MTX_LOGI("Radical length is greater than %d\n",soclen);
	if (!Head)
	    return 0;	
     }

    //----------------------------------------------------------------------------

    if (soc2 == NULL)
    {
	MTX_LOGI("Radical length is smaller than %d, there are no vectors "
		"in the %dth Head\n",Len,Len);
	return 0;
    }

   // computing the <Len>th radical

    rad2 = matNullSpace__(matTransposed(soc2));	// rad^Len given in rad^(Len - 1)

    Dualize(Rep);	// the action on rad^(Len - 1)

    // calculates the vecors generating the irreducibles lying in rad^Len
    rep = wgAlloc(Rep);
    for(j = 0; j < Info->nCf; j++)
    {
	Matrix_t *word = wgMakeWord(rep,Info->Cf[j].peakWord);
	Matrix_t *seed2;

	// makes the iterated nullspace of the peakword
	matInsert_(word, Info->Cf[j].peakPol);
	if ((seed = matAlloc(ffOrder, 0, 0)) == NULL)
	    return 1;
	for (seed2 = matNullSpace(word); seed->nor < seed2->nor; 
	    seed2 = matNullSpace(matMul(word, word)))
	{
	    matFree(seed);
	    seed = seed2;
	}
	matFree(seed2);
	matFree(word);

	if (seed->nor > 0)
	{
	    Matrix_t *sec = intersect(seed, rad2); // the nullspace in rad^2
	    matPivotize(sec);
	    matClean(seed,sec);
	    matFree(sec);
	    if (emb != NULL)
		matMul(seed, emb); // embedding into the original module
	    mat = seed;
	}
	else
	{
	    int noc = (emb == NULL) ? Rep->Gen[0]->noc : emb->noc;
	    matFree(seed);
	    if ((mat = matAlloc(ffOrder, 0, noc)) == NULL)
	    {
		mtxAbort(MTX_HERE,"Cannot allocate matrix");
		return 1;
	    }
	}

	sprintf(name,"%s%s.h%d", Name, latCfName(Info,j), Len);
	matSave(mat,name);
/*	matPrint(name, mat);*/
	matFree(mat);
    }
    wgFree(rep);
    latDestroy(Info);
    return(0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

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
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
