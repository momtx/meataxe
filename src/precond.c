////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Precondensation of tensor products
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

char TkiName[LAT_MAXBASENAME];	/* Base name for .tki file */
TkData_t TKInfo;			/* Data from .tki file */
Lat_Info InfoM, InfoN;		/* Data from .cfinfo files */
int CfIsLinked[LAT_MAXCF];	/* Constituent (of N) is already linked */
Matrix_t *Trans[LAT_MAXCF];	/* */
int opt_s = 0;			/* Do NOT use standard generators */


static MtxApplicationInfo_t AppInfo = { 
"precond", "Precondensation",
"SYNTAX"
"    precond " MTX_COMMON_OPTIONS_SYNTAX " <Info> <M> <N>\n"
"\n"
"ARGUMENTS\n"
"    <Info> .................. Name for condensation info file\n"
"    <M> ..................... First module (left factor), semisimple\n"
"    <N> ..................... Second module (right factor), semisimple\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"\n"
"FILES\n"
"    <M>.cfinfo .............. I  Info file for M (produced by CHOP/PWKOND)\n"
"    <N>.cfinfo .............. I  Info file for N (produced by CHOP/PWKOND)\n"
"    <M/N><Cf>.std.{1,2...} .. I  Standard generators for each constituent\n"
"    <Info>.tki .............. O  Tensor condensation info file\n"
"    <Info>.q.{1,2...} ....... O  Embeddings for each constituent\n"
"    <Info>.p.{1,2...} ....... O  Projections for each constituent\n"
};

static MtxApplication_t *App = NULL;



/* --------------------------------------------------------------------------
   Init() - Initialize everything

   Description:
     This function initializes global variables, processes command line
     options and arguments, and reads the .cfinfo and .tki files.

   Arguments:
     <argc>: Number of arguments
     <argv>: Command line arguments
   -------------------------------------------------------------------------- */

static int Init(int argc, const char **argv)

{
    int i;
    int fl1, nor1, noc1, fl2, nor2, noc2;
    char fn[LAT_MAXBASENAME+20];

    /* Initialize global data
       ---------------------- */
    memset(&TKInfo,0,sizeof(TKInfo));
    memset(TkiName,0,sizeof(TkiName));
    memset(CfIsLinked,0,sizeof(CfIsLinked));

    /* Initialize the MeatAxe library
       ------------------------------ */
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    opt_s = AppGetOption(App,"-s");
    if (AppGetArguments(App,3,3) != 3)
	return -1;
    strcpy(TkiName,App->ArgV[0]);

    /* Read info files
       --------------- */
    strncpy(TKInfo.NameM,App->ArgV[1],sizeof(TKInfo.NameM));
    strncpy(TKInfo.NameN,App->ArgV[2],sizeof(TKInfo.NameN));
    if (Lat_ReadInfo(&InfoM,App->ArgV[1]) != 0)
    {
	MTX_ERROR1("Error reading %s.cfinfo",App->ArgV[1]);
	return -1;
    }
    if (Lat_ReadInfo(&InfoN,App->ArgV[2]) != 0)
    {
	MTX_ERROR1("Error reading %s.cfinfo",App->ArgV[2]);
	return -1;
    }


    /* Initialize the TKInfo structure
       ------------------------------- */
    TKInfo.NCf = 0;
    for (i = 0; i < LAT_MAXCF; ++i)
	TKInfo.CfIndex[0][i] = TKInfo.CfIndex[1][i] = -1;

    /* Some additional checks on input files
       ------------------------------------- */
    if (InfoM.Field != InfoN.Field)
    {
	MTX_ERROR4("Incompatible representations: %s is over GF(%d), "
	   "%s is over GF(%d)",InfoM.BaseName,InfoM.Field,InfoN.BaseName,
	   InfoN.Field);
	return -1;
    }
    if (InfoM.NGen != InfoN.NGen)
    {
	MTX_ERROR4("Incompatible representations: %s has %d generators, "
	   "%s has %d generators",InfoM.BaseName,InfoM.NGen,InfoN.BaseName,
	   InfoN.NGen);
	return -1;
    }

    /* Print start message
       ------------------- */
    sprintf(fn,"%s.1",InfoM.BaseName);
    fclose(FfReadHeader(fn,&fl1,&nor1,&noc1));
    sprintf(fn,"%s.1",InfoN.BaseName);
    fclose(FfReadHeader(fn,&fl2,&nor2,&noc2));
    MESSAGE(0,("Beginning pre-condensation of dimension %d x %d = %d\n",
	nor1,nor2,nor1*nor2));
    return 0;
}



/* --------------------------------------------------------------------------
   IsDual() - Find matching constituent

   Description:
     This function decides if a given contituent of M is dual to a given 
     contituent of N.

   Arguments:
     <mj>: Index of constituent of M.
     <rep_m>: Generators of the constituent of M.
     <nj>: Index of constituent of N.

   Return:
     1 = the constituents are dual, 0 = they are not dual.
   -------------------------------------------------------------------------- */

static int IsDual(int mj, MatRep_t *rep_m, int nj)

{
    MatRep_t *rep_n;	    /* Generators of constituent of N */
    int result;
    CfInfo *minfo = InfoM.Cf + mj;

    /* First check: Dimensions and splitting field must match
       ------------------------------------------------------ */
    if (InfoN.Cf[nj].dim != minfo->dim || InfoN.Cf[nj].spl != minfo->spl)
	return 0;

    /* Read the (contragrediate) generators and compare
       ------------------------------------------------ */
    MESSAGE(2,(" (%s%s)",InfoN.BaseName,Lat_CfName(&InfoN,nj)));
    rep_n = Lat_ReadCfGens(&InfoN,nj,LAT_RG_INVERT|LAT_RG_TRANSPOSE
	| (InfoN.Cf[nj].peakword >= 0 ? LAT_RG_STD : 0));
    
    result = IsIsomorphic(rep_m,minfo,rep_n,Trans + TKInfo.NCf,
	minfo->peakword >= 0);
    MrFree(rep_n);
    return result;
}




/* --------------------------------------------------------------------------
   FindConstituentInN() - Find matching constituent

   Description:
     This function finds, for a given constituent of M, the corresponding
     constituent in N (if any exists). The generators in <mgen> are supposed
     to be in standard basis with respect to the constituent's idword.

   Arguments:
     <mj>: Index of the constituent of M.
     <rep_m>: Generators of the constituent of M.

   Return:
     Index of corresponding constituent of N, or -1 if not found.
   -------------------------------------------------------------------------- */

static int FindConstituentInN(int mj, MatRep_t *rep_m)

{
    int nj;

    /* Main loop: for each constituent of N
       ------------------------------------ */
    for (nj = 0; nj < InfoN.NCf; ++nj)
    {
	/* Discard constituents which are already linked
	   --------------------------------------------- */
	if (CfIsLinked[nj])
	   continue;

	/* Compare with <nj>-th constituent of N
	   ------------------------------------- */
	if (IsDual(mj,rep_m,nj))
	{
	    CfIsLinked[nj] = 1;	    /* Mark as linked */
	    return nj;
	}
    }

    /* Corresponding constituent was not found
       --------------------------------------- */
    return -1;
}



/* --------------------------------------------------------------------------
   MkEndo() - Calculate endomorphism ring

   Description:
     This function calculates, for a given irreducible constituent, a
     basis for the endomoprhism ring.

   Arguments:
     <rep>: Generators.
     <cf>: Pointer to constituent information.
     <endo>: Pointer to an array where the basis of End(V) is stored.
     <maxendo>: Maximum number of matrices that can be stored in <endo>.
     <basis>: The 'E-basis' is stored here
   -------------------------------------------------------------------------- */

static void MkEndo(const MatRep_t *rep, const CfInfo *cf,
    Matrix_t **endo, int maxendo)

{
    Matrix_t *pw, *nsp;
    WgData_t *wg;

    MTX_VERIFY(maxendo >= cf->spl);

    /* Make the peak word kernel
       ------------------------- */
    wg = WgAlloc(rep);
    pw = WgMakeWord(wg,cf->idword);
    WgFree(wg);
    nsp = MatNullSpace__(MatInsert(pw,cf->idpol));
    MTX_ASSERT(nsp->Nor == cf->spl);
    MatFree(pw);

    /* Calculate a basis of the the endomorphism ring
       ---------------------------------------------- */
    const int i = MakeEndomorphisms(rep,nsp,endo);
    (void)i;
    MTX_ASSERT(i == 0);

    MatFree(nsp);
}




/* --------------------------------------------------------------------------
   MakeQ() - Calculate the Q matrix

   Description:
     This function calculates the embedding of the condensed tensor product
     (VxW)e into the tensor product VxW for one irreducible constituent.

   Arguments:
     <n>: Constituent index.
     <spl>: Splitting field degree.
     <endo>: Basis of endomorphism ring.
   -------------------------------------------------------------------------- */

static void MakeQ(int n, int spl, const Matrix_t **endo)

{		  
    int i;
    int dim = endo[0]->Nor;
    Matrix_t *q = MatAlloc(endo[0]->Field,spl,dim*dim);
    char fn[200];
    for (i = 0; i < spl; ++i)
    {
	int j;
	Matrix_t *y = MatInverse(Trans[n]), *x;
	MatMul(y,endo[i]);
	x = MatTransposed(y);
	MatFree(y);
	for (j = 0; j < dim; ++j)
	    MatCopyRegion(q,i,j * dim,x,j,0,1,-1);
	MatFree(x);
    }
    sprintf(fn,"%s.q.%d",TkiName,n+1);
    MESSAGE(2,("Writing %s\n",fn));

#if 0
    if (InfoM.Cf[TKInfo.CfIndex[0][n]].peakword < 0)
	MatMulScalar(q,FF_ZERO);
#endif

    MatSave(q,fn);
    MatFree(q);
}



static FEL MatProd(Matrix_t *a, Matrix_t *b)

{
    FEL f = FF_ZERO;
    int i;

    FfSetNoc(a->Noc);
    for (i = 0; i < a->Nor; ++i)
	f = FfAdd(f,FfScalarProduct(MatGetPtr(a,i),MatGetPtr(b,i)));
    return f;
}





/* --------------------------------------------------------------------------
   MakePQ() -

   Description:

   Arguments:
     <n>: Constituent index
     <mj>: Index of corresponding constituent of M.
     <nj>: Index of corresponding constituent of N.
   -------------------------------------------------------------------------- */

#define MAXENDO 20

static void MakePQ(int n, int mj, int nj)

{
    MatRep_t *rep_m;
    Matrix_t *estar[MAXENDO], *endo[MAXENDO], *e, *ei;
    char fn[200];
    int dim = InfoM.Cf[mj].dim;
    int spl = InfoM.Cf[mj].spl;
    int i;
    Matrix_t *p;

    MESSAGE(1,("Condensing %s%s x ",InfoM.BaseName,Lat_CfName(&InfoM,mj)));
    MESSAGE(1,("%s%s, [E:k]=%d\n",InfoN.BaseName,Lat_CfName(&InfoN,nj),spl));

    /* Read the generators for the constituent of M and make the
       endomorphism ring.
       --------------------------------------------------------- */
    rep_m = Lat_ReadCfGens(&InfoM,mj,InfoM.Cf[mj].peakword >= 0 ? LAT_RG_STD : 0);
    MESSAGE(2,("Calculating endomorphism ring\n"));
    MkEndo(rep_m,InfoM.Cf + mj,endo,MAXENDO);
    MrFree(rep_m);

    /* Calculate the Q matrix
       ---------------------- */
    MESSAGE(2,("Calculating embedding of E\n"));
    MakeQ(n,spl,(const Matrix_t **)endo);
    
    /* Calculate the E* matrices
       Note: We should use the symmetry under i<-->k here!
       --------------------------------------------------- */
    MESSAGE(2,("Calculating projection on E\n"));
    MESSAGE(2,("   E* matrices\n"));
    e = MatAlloc(FfOrder,spl,spl);
    for (i = 0; i < spl; ++i)
    {
	PTR pptr = MatGetPtr(e,i);
	int k;
	for (k = 0; k < spl; ++k)
	{
	    FEL f;
	    Matrix_t *x = MatDup(endo[i]);  
	    MatMul(x,endo[k]);
	    f = MatTrace(x);
	    FfInsert(pptr,k,f);
	    MatFree(x);
	}
    }
    ei = MatInverse(e);
    MatFree(e);

    for (i = 0; i < spl; ++i)
    {
	int k;
	PTR p;
	estar[i] = MatAlloc(FfOrder,dim,dim);
	p = MatGetPtr(ei,i);
	for (k = 0; k < spl; ++k)
	    MatAddMul(estar[i],endo[k],FfExtract(p,k));
    }
    MatFree(ei);

    /* Transpose the E* matrices. This simplifies the 
       calculation of tr(z E*) below.
       ----------------------------------------------- */
    MESSAGE(2,("   Transposing E* matrices\n"));
    for (i = 0; i < spl; ++i)
    {
	Matrix_t *x = MatTransposed(estar[i]);
	MatFree(estar[i]);
	estar[i] = x;
    }

    /* Calculate the P matrix
       ---------------------- */
    MESSAGE(2,("   P matrix\n"));
    p = MatAlloc(FfOrder,dim*dim,spl);
    for (i = 0; i < dim; ++i)
    {
	int j;
	for (j = 0; j < dim; ++j)
	{
	    int r;
	    PTR pptr = MatGetPtr(p,i*dim + j);
	    Matrix_t *x = MatAlloc(FfOrder,dim,dim);
	    MatCopyRegion(x,0,i,Trans[n],0,j,dim,1);
	    for (r = 0; r < spl; ++r)
	    {
		FEL f = MatProd(x,estar[r]);
		FfInsert(pptr,r,f);
	    }
	    MatFree(x);
	}
    }

    sprintf(fn,"%s.p.%d",TkiName,n+1);
    MESSAGE(2,("Writing %s\n",fn));
#if 0
    if (InfoM.Cf[mj].peakword < 0)
	MatMulScalar(p,FF_ZERO);
#endif
    MatSave(p,fn);

    /* Clean up
       -------- */
    MatFree(p);
    for (i = 0; i < spl; ++i)
        MatFree(endo[i]);
}




/* --------------------------------------------------------------------------
   CalcDim() - Calculate dimension of condensed module

   Description:
     This function calculates the dimension of the condensed tensor product.

   Remarks:
     The dimension is stored into <TKInfo.Dim>.
   -------------------------------------------------------------------------- */

static void CalcDim()

{
    int i;

    TKInfo.Dim = 0;
    for (i = 0; i < TKInfo.NCf; ++i)
    {
	int m = TKInfo.CfIndex[0][i];
	int n = TKInfo.CfIndex[1][i];
	TKInfo.Dim += InfoM.Cf[m].mult * InfoN.Cf[n].mult * InfoM.Cf[m].spl;
    }
    MESSAGE(0,("Dimension of condensed module = %d\n",TKInfo.Dim));
}


/* --------------------------------------------------------------------------
   main() - Program entry point
   -------------------------------------------------------------------------- */

int main(int argc, const char **argv)

{
    int mj;

    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return -1;
    }

    /* Main loop: for all constituents of M
       ------------------------------------ */
    for (mj = 0; mj < InfoM.NCf; mj++)
    {
	MatRep_t *rep_m;	    /* Generators for const. of M */
	int nj;

	if (InfoM.Cf[mj].peakword < 0)
	{
	    MESSAGE(0,("WARNING: No peak word word available for %s%s\n",
		InfoM.BaseName,Lat_CfName(&InfoM,mj)));
	}
	MESSAGE(0,("%s%s ",InfoM.BaseName,Lat_CfName(&InfoM,mj)));

	/* Read the generators for the <mj>-th contituent of M, and find
	   the corresponding (=contragredient) constituent in N.
	   ------------------------------------------------------------- */
	rep_m = Lat_ReadCfGens(&InfoM,mj,InfoM.Cf[mj].peakword >= 0 ? LAT_RG_STD : 0);
	nj = FindConstituentInN(mj,rep_m);

	/* Calculate the P and Q matrix for this constituent
	   ------------------------------------------------- */
	if (nj >= 0)
	{
	    MESSAGE(0,(" <--> %s%s\n",InfoN.BaseName,
		Lat_CfName(&InfoN,nj)));
	    TKInfo.CfIndex[0][TKInfo.NCf] = mj;
	    TKInfo.CfIndex[1][TKInfo.NCf] = nj;
	    MakePQ(TKInfo.NCf,mj,nj);
	    TKInfo.NCf++;
	}
	else
	    MESSAGE(0,(" not found in %s\n",InfoN.BaseName));

	/* Clean up
	   -------- */
	MrFree(rep_m);
    }

    CalcDim();				/* Calculate dimension */
    TK_WriteInfo(&TKInfo,TkiName);	/* Write .tki file */
    if (App != NULL) AppFree(App);
    return 0;
}


/**
@page prog_precond precond - Precondensation of Tensor Products"

@section precond_syntax Command Line
<pre>
precond @em Options @em Info @em M @em N
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par @em Info
  Name of the tensor condensation data file.
@par @em M
  Name of first module (left factor).
@par @em N
  Name of second module (right factor).

@section precond_inp Input Files
@par @em M.cfinfo, @em N.cfinfo
  Constituent information.
@par @em MCf.std.1, @em MCf.std.2, ..., @em NCf.std.1, @em NCf.std.2, ...
  Standard generators of the condensation subgroup H for each constituent.

@section precond_out Output Files
@par @em Info.tki
  Tensor condensation info file.
@par @em Info.q.1, @em Info.q.2, ...
  Embeddings for each constituent.
@par @em Info.p.1, @em Info.p.2, ...
  Projections for each constituent.


@section precond_desc Description
This program performs two tasks:
- It compares the irreducible constituents of M<sub>H</sub> and N<sub>H</sub>,
  and finds all pairs (S<sub>i</sub>,T<sub>j</sub>) of constituents where
  S<sub>i</sub>≅T<sub>j</sub><sup>⋆</sup>.
- For each pair (S,T) of constituents found in step 1, the program
  calculates the embedding of (S⊗T)e<sub>H</sub> into S⊗T as an 
  direct summand, and the corresponding projection of S⊗T onto
  (S⊗T)e<sub>H</sub>.
If there is no peak word for a constituent, @b precond will issue a warning
but continue. However, the P and Q matrices for this constituent are zero.

@section precond_impl Implementation Details
Step 1, matching of constituents, is implemented in the same way as in @ref prog_chop "chop"
and @ref prog_cfcomp "cfcomp", i.e., by using the standard basis with respect to identifying 
words.
Step 2 is based on two observations:
- (A):
    V⊗V<sup>*</sup>≅Hom<sub>k</sub>(V,V)$, and (S⊗T)e_H≅End<sub>kH</sub>(V)
    as $kH$-Modules.
- (B):
    There is a natural, H-invariant non-degenerate scalar product on 
    Hom<sub>k</sub>(V,V), given by Γ(φ,ψ)=Trace(φ∘ψ).

From (A) it is clear that calculating the embedding of
(S⊗T)e<sub>H</sub> into S⊗T is equivalent to computing a basis
of End<sub>kH</sub>(V).
The latter is easily accomplished using the peak word of V.
As a consequence of the second observation, there is a natural one-to-one 
correspondence between H-invariant linear forms on Hom_k(V,V)
and End<sub>kH</sub>(V), which is used to calculate the projection from
Hom<sub>k</sub>(V,V) on End<sub>kH</sub>(V).

More details on the algorithm used in Step 2 can be found in @ref Ri98 "[Ri98]".

**/
