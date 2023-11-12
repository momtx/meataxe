////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Condensation of tensor products.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

/* --------------------------------------------------------------------------
   Global data
   -------------------------------------------------------------------------- */



static const char *TkiName;		/* Base name for .tki file */
static const char *ResultName;		/* Base name for result */
static TkData_t TKInfo;			/* Data from .tki file */
static Lat_Info InfoM, InfoN;	        /* Data from .cfinfo files */
static const char *AName;		/* Left factor */
static const char *BName;	        /* Right factor */
static int NGen = 2;			/* No. of generators for A and B */
static Matrix_t *SsBasisM, *SsBasisN;	/* Semisimplicity basis */
static Matrix_t *SsBasisMi, *SsBasisNi;
static Matrix_t *Q[LAT_MAXCF];		/* Q matrices (embeddings) */
static Matrix_t *P[LAT_MAXCF];		/* P matrices (projections) */
static int WriteGenerators = 0;		/* -t: Write transformed gens of A, B */
static int NoBasisChange = 0;		/* -n: No basis change */

static MtxApplicationInfo_t AppInfo = { 
"tcond", "Condense tensor product", 
"\n"
"SYNTAX\n"
"    tcond [-QVt] [-T <MaxTime>] [-g <NGen>] <Info> <A> <B> <Result>\n"
"\n"
"ARGUMENTS\n"
"    <Info> .................. Base name of .tki file\n"
"    <A>, <B> ................ Representations of G\n"
"    <Result> ................ Name for condensed representation\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -g <NGen> ............... Set number of generators (default: 2)\n"
"    -t ...................... Write transformed generators of <A> and <B>\n"
"    -n ...................... No basis change. Assume that generators on <A>\n"
"                              and <B> are already in semisimplicity basis.\n"
"\n"
"FILES\n"
"    <Info>.tki .............. I Tensor condensation info, made by PRECOND\n"
"    <A>.{1,2...} ............ I Generators of left factor\n"
"    <B>.{1,2...} ............ I Generators of right factor\n"
"    <A>.ssb ................. I Semisimplicity basis for A, made by PWKOND\n"
"    <B>.ssb ................. I Semisimplicity basis for B, made by PWKOND\n"
"    <Info>.q.{1,2...} ....... I Basis matrices for constituents\n"
"    <Info>.p.{1,2...} ....... I Projection matrices for constituents\n"
"    <Result>.{1,2...} ....... O Condensed matrices\n"
"    <A>.ss.{1,2...} ......... O Transformed generators (with -t)\n"
"    <B>.ss.{1,2...} ......... O Transformed generators (with -t)\n"
};
static MtxApplication_t *App = NULL;




static void MakeInvertible(Matrix_t *mat, const char *fn)
{
    Matrix_t *dup = matDup(mat);
    PTR x;
    int i, k;

    matEchelonize(dup);
    k = dup->nor;
    if (k < mat->nor)
    {
	MESSAGE(0,("WARNING: %s: %d basis vectors are missing, "
	    "using random vectors\n",fn,mat->nor - k));
    }
    for (i = 0, x = mat->data; i < mat->nor; ++i, ffStepPtr(&x, mat->noc))
    {
	FEL f;
	if (ffFindPivot(x,&f, mat->noc) == MTX_NVAL)
	{
	    ffInsert(x,dup->pivotTable[k],FF_ONE);
	    ++k;
	}
    }
    MTX_ASSERT(k == mat->nor);
    matFree(dup);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void Init(int argc, char** argv)
{
   char fn[200];
   int i;

   // Process command line options
   App = appAlloc(&AppInfo, argc, argv);
   NGen = appGetIntOption(App, "-g", 2, 1, 100);
   WriteGenerators = appGetOption(App, "-t --write-generators");
   NoBasisChange = appGetOption(App, "-n --no-basis-change");
   if (WriteGenerators && NoBasisChange) {
      mtxAbort(NULL, "'-t' and '-n' cannot be used together");
   }

   // Process command line arguments
   appGetArguments(App, 4, 4);
   TkiName = App->argV[0];
   AName = App->argV[1];
   BName = App->argV[2];
   ResultName = App->argV[3];

   // Read .tki file
   {
      int ctx = mtxBegin(MTX_HERE, "HINT: did you run 'precond'?");
      tkReadInfo(&TKInfo, TkiName);
      if (TKInfo.dim <= 0) {
         mtxAbort(MTX_HERE, "No dimension found in .tki file");
      }
      mtxEnd(ctx);
   }

   latReadInfo(&InfoM, TKInfo.nameM);
   latReadInfo(&InfoN, TKInfo.nameN);

   // Some checks on info file data
   if (InfoM.field != InfoN.field) {
      mtxAbort(MTX_HERE, "Different fields in .cfinfo files");
   }
   if (InfoN.NGen != InfoM.NGen) {
      mtxAbort(MTX_HERE, "Different number of generators in %s and %s",
         InfoM.BaseName, InfoN.BaseName);
   }

   // Read the semisimplicity bases.
   if (!NoBasisChange) {
      int ctx = mtxBegin(MTX_HERE, "HINT: did you run 'pwkond -tb'?");
      MESSAGE(1, ("Reading and inverting semisimplicity bases\n"));
      sprintf(fn, "%s.ssb", TKInfo.nameM);
      SsBasisM = matLoad(fn);
      MakeInvertible(SsBasisM, fn);
      SsBasisMi = matInverse(SsBasisM);
      if (strcmp(AName, BName)) {
         sprintf(fn, "%s.ssb", TKInfo.nameN);
         SsBasisN = matLoad(fn);
         MakeInvertible(SsBasisN, fn);
         SsBasisNi = matInverse(SsBasisN);
      }
      else {
         SsBasisN = SsBasisM;
         SsBasisNi = SsBasisMi;
      }
      mtxEnd(ctx);
   }

   // Read P and Q matrices.
   for (i = 0; i < TKInfo.nCf; ++i) {
      int ctx = mtxBegin(MTX_HERE, "HINT: did you run 'precond'?");
      int spl = InfoM.Cf[TKInfo.cfIndex[0][i]].spl;
      int tdim = InfoM.Cf[TKInfo.cfIndex[0][i]].dim;
      int f;
      tdim *= tdim;
      sprintf(fn, "%s.p.%d", TkiName, i + 1);
      P[i] = matLoad(fn);
      f = ffOrder;
      if (P[i]->field != f || P[i]->noc != spl || P[i]->nor != tdim) {
         mtxAbort(MTX_HERE, "%s: Invalid P matrix", fn);
      }
      sprintf(fn, "%s.q.%d", TkiName, i + 1);
      Q[i] = matLoad(fn);
      if (Q[i]->field != f || Q[i]->nor != spl || Q[i]->noc != tdim) {
         mtxAbort(MTX_HERE, "%s: Invalid Q matrix", fn);
      }
      mtxEnd(ctx);
   }
}



/* --------------------------------------------------------------------------
   FirstRow() - Get first row in Q matrix

   Description:
     This function calculates the index of the first row in the Q matrix 
     which belongs to a given constituent. The constituent is identified by
     it isomorphism type <cf>, and, if the constituent occurs more than once,
     an additional index, <k>, running from 0 to $m-1$, where $m$ is the 
     multiplicity.
     All Indexes are 0-based. 

   Arguments:
     <info>: Pointer to constituent information
     <cf>: Constituent index
     <k>: Counter for multiple copies of the same constituent.

   Return:
     Index of first row belonging to the constituent.
   -------------------------------------------------------------------------- */

static int FirstRow(Lat_Info *info, int cf, int k)

{
    int ind = 0;
    int i;

    MTX_ASSERT(cf >= 0 && cf < info->nCf);
    MTX_ASSERT(k >= 0 && k < info->Cf[cf].mult);

    /* Constituents before <cf> consume dim*mult rows
       ---------------------------------------------- */
    for (i = 0; i < cf; ++i)
	ind += info->Cf[i].dim * info->Cf[i].mult;

    /* Copies of constituent <cf> before <k> consume dim rows
       ------------------------------------------------------ */
    ind += info->Cf[cf].dim * k;

    /* Return 0-based row index
       ------------------------ */
    return ind;
}








/*----------------------------------------------------------------------------*/

static void gemap(Matrix_t *conma, Matrix_t *q, Matrix_t *mrow, Matrix_t *nrow)

  /* the vectors from "q" are mapped under g*e */

{
    int j;
    int bcol = 0;


    /* For each irreducible constituent I
       ---------------------------------- */
    for (j = 0; j < TKInfo.nCf; ++j)
    {
	int cfm = TKInfo.cfIndex[0][j];		/* Index of constituent in M */
	int cfn = TKInfo.cfIndex[1][j];		/* Index of constituent in M */
	int d = InfoM.Cf[cfm].dim;		/* Dimension */
	int mj;
	
	/* For each copy of I in M
	   ----------------------- */
        for (mj = 0; mj < InfoM.Cf[cfm].mult; ++mj)
        {
	    long mstart = FirstRow(&InfoM,cfm,mj);
	    Matrix_t *mop = matCut(mrow,0,mstart,-1,d);
	    int nj;

	    /* For each copy of I in N
	       ----------------------- */
            for (nj = 0; nj < InfoN.Cf[cfn].mult; nj++)
            {   
		long nstart = FirstRow(&InfoN,cfn,nj);
                Matrix_t *nop = matCut(nrow,0,nstart,-1,d);
		Matrix_t *image = TensorMap(q,mop,nop);
                matFree(nop);
                matMul(image, P[j]);   /* projection */
		matCopyRegion(conma,0,bcol,image,0,0,-1,-1);
                matFree(image);
                bcol += P[j]->noc;
            }
	    matFree(mop);
        }
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Condense one generator

static void condenseMat(int gen)
{
    char resname[LAT_MAXBASENAME + 20];	/* Output file name */
    char aname[LAT_MAXBASENAME + 20];	/* Generator on M */
    char bname[LAT_MAXBASENAME + 20];	/* Generator on M */
    int cf;				/* Constituent index */
    Matrix_t *mmat, *nmat, *x;		/* The generator on M, N */


    /* Make file names 
       --------------- */
    sprintf(resname,"%s.%d",ResultName,gen+1);
    sprintf(aname,"%s.%d",AName,gen+1);
    sprintf(bname,"%s.%d",BName,gen+1);
    MESSAGE(0,("Condensing %s x %s --> %s\n",aname,bname,resname));

    /* Load the generator on M and N
       ----------------------------- */
    mmat = matLoad(aname);
    if (strcmp(aname,bname))
	nmat = matLoad(bname);
    else
	nmat = mmat;

    /* Change to semisimplicity basis.
       ------------------------------- */
    if (!NoBasisChange)
    {
        MESSAGE(1,("  Changing basis\n"));
        x = matDup(SsBasisM);
        matMul(x,mmat);
        
        matMul(x,SsBasisMi);
        matFree(mmat);
        mmat = x;

        if (strcmp(aname,bname))
        {
	    x = matDup(SsBasisN);
	    matMul(x,nmat);
	    matMul(x,SsBasisNi);
            matFree(nmat);
	    nmat = x;
        }
        else
	    nmat = mmat;
    }

    if (WriteGenerators)
    {
	char fn[200];
        sprintf(fn,"%s.ss.%d",AName,gen+1);
	matSave(mmat,fn);
	if (strcmp(aname,bname))
        {
	    sprintf(fn,"%s.ss.%d",BName,gen+1);
	    matSave(nmat,fn);
	}
    }

    /* Open the output file
       -------------------- */
    MESSAGE(1,("  Beginning condensation\n"));
    MtxFile_t* resultFile = mfCreate(resname,ffOrder,TKInfo.dim,TKInfo.dim);

    /* Main loop: for each constituent
       ------------------------------- */
    for (cf = 0; cf < TKInfo.nCf; ++cf)
    {
	int cfm = TKInfo.cfIndex[0][cf];    /* Index in M */
	int cfn = TKInfo.cfIndex[1][cf];    /* Index in N */
	int rownb = InfoM.Cf[cfm].dim;	/* Number of rows to extract */
	int mi;				/* Counter for copies of this const. */

        MESSAGE(2,("  Processing %s",latCfName(&InfoM,cfm)));
        MESSAGE(2,(" x %s\n",latCfName(&InfoN,cfn)));

        for (mi = 0; mi < InfoM.Cf[cfm].mult; ++mi)
        {
	    Matrix_t *mrow;
	    int ni;
	    int firstrow = FirstRow(&InfoM,cfm,mi);
	    mrow = matCutRows(mmat,firstrow,rownb);

	    MESSAGE(3,("  "));
            for (ni = 0; ni < InfoN.Cf[cfn].mult; ++ni)
            {
		Matrix_t *nrow, *condmat;

		firstrow = FirstRow(&InfoN,cfn,ni);
		nrow = matCutRows(nmat,firstrow,rownb);
		condmat = matAlloc(ffOrder,Q[cf]->nor,TKInfo.dim);
		if (condmat == NULL)
		{
		    mtxAbort(MTX_HERE,"Cannot allocate %dx%d matrix",
			Q[cf]->nor,TKInfo.dim);
		}

	    	MESSAGE(3,(" %dx%d",mi,ni));
                gemap(condmat,Q[cf],mrow,nrow);
                
                /* write result */
        	mfWriteRows(resultFile,condmat->data,condmat->nor,condmat->noc);
                matFree(condmat);
                matFree(nrow);
            }
	    MESSAGE(3,("\n"));
            matFree(mrow);
        }   
    }   

    mfClose(resultFile);

    /* Clean up
       -------- */
    matFree(mmat);
    if (strcmp(aname,bname))
	matFree(nmat);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    Init(argc,argv);
    for (int i = 0; i < NGen; ++i)
    {
	// TODO: context "Condensation of %s.%d x %s.%d",AName,i+1,BName,i+1);
	condenseMat(i);
    }
    appFree(App);
    return 0;
}  



/**
@page prog_tcond tcond - Tensor Product Condensation

@section tcond_syntax Command Line
<pre>
tcond [@em Options] [-nt] [-T @em MaxTime] [-g @em NGen] @em Info @em M @em N @em Result
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -g @em NGen
  Set the number of generators. Default is 2.
@par -n
  No basis change. Assume that the generators on M and N are already given with
  respect to the semisimplicity basis.
@par -t
  Write out transformed generators for M and N.
@par @em Info
  Information produced by @ref prog_precond "precond".
@par @em M
  Name of the first representation (left factor).
@par @em N
  Name of the second representation (right factor).
@par @em Result
  Condensed representation.

@section tcond_inp Input Files
@par @em Info.tki
  Tensor condensation info, made by @ref prog_precond "precond".
@par @em M.1, @em M.2, ...
  Generators of left factor.
@par @em N.1, @em N.2, ...
  Generators of right factor.
@par @em M.ssb
  Semisimplicity basis for M, made by @ref prog_pwkond "pwkond".
@par @em N.ssb
  Semisimplicity basis for N, made by @ref prog_pwkond "pwkond".
@par @em Info.q.1, @em Info.q.2, ...
  Basis matrices for constituents.
@par @em Info.p.1, @em Info.p.2, 
  Projection matrices for constituents.

@section tcond_out Output Files
@par @em Result.1, @em Result.2, ...
  Condensed matrices.
@par @em M.ss.1, @em M.ss.2, ...
  Transformed generators (with -t)
@par @em N.ss.1, @em N.ss.2, ...
  Transformed generators (with -t)

@section tcond_desc Description
This program performs the final steps of the tensor condensation procedure.
It calculates, for one or more elements a₁,a₂,…∊A,
the action of e<sub>H</sub>a<sub>i</sub>e<sub>H</sub> on the condensed tensor
product (M⊗N)e<sub>H</sub>.


As input, the program expects the action of a<sub>i</sub> on M and N with
respect to the same basis as the generators of the condensation 
subgroup H that were fed into @ref prog_precond "precond" before.
The program also needs the semisimplicity basis calculated by
@ref prog_pwkond "pwkond", and the P and Q matrices calculated by @ref prog_precond "precond".

If the generators are already given with repect to the semisimplicity
basis, you can use the "-n" option to tell @b tcond to skip the basis
change.

The output are @em NGen matrices describing the action of
e_<sub>H</sub>a<sub>i</sub>e<sub>H</sub> on
(M⊗N)e<sub>H</sub>. These matrices are written to @em Result.1, @em Result.2 ...
If you use the "-t" option, @b tcond also calculates 
the action of a<sub>i</sub> on M and N with respect to the semisimplicity
basis. This option cannot be used together with "-n".

The following sequence of commands shows the complete procedure for
condensing a tensor product. To make things more simple, we assume that
M=N. The condensation subgroup shall be given by three generators in the
files "sub.1", "sub.2", and "sub.3".
The generators of the group shall be "g.1" and "g.2".
<pre>
chop -g 3 sub
pwkond -tb sub
precond tp sub sub
tcond -g 2 tp g g result
</pre>
After these commands are completed, the action of the condensed generators
is in "result.1", "result.2", and "result.3".

@section tcond_impl Implementatin Details
The algorithm used by this program is described in @ref Wie94 "[Wie94]".


*/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
