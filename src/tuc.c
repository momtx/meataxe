////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - This program uncondenses one or more vectors in (M x N)e.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

/* --------------------------------------------------------------------------
   Global data
   -------------------------------------------------------------------------- */


const char *TkiName;		    // Base name for .tki file
TkData_t TkInfo;		    // Data from .tki file
LatInfo_t* InfoM;		    // Data from .cfinfo files
LatInfo_t* InfoN;		    // Data from .cfinfo files
Matrix_t *CondMat;		    // Condensed vectors (read from file)
Matrix_t *UncondMat;		    // Uncondensed vectors (calculated here)
const char *UncondName;		    // Name of output file
long DimM, DimN;		    // Dimensions of M and N
Matrix_t *QMat[LAT_MAXCF];	    // Q matrices

static MtxApplicationInfo_t AppInfo = { 
"tuc", "Tensor Uncondense", 
"\n"
"SYNTAX\n"
"    tuc " MTX_COMMON_OPTIONS_SYNTAX " <info> <cond> <uncond>\n"
"\n"
"ARGUMENTS\n"
"    <info> .................. Tensor condensation file name\n"
"    <cond> .................. Vectors to uncondense\n"
"    <uncond>................. Condensed vectors\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"\n"
"FILES\n"
"    <info>.tki .............. I Tensor condense info file\n"
"    <M>.cfinfo .............. I Constituent info file for left factor\n"
"                                (name taken from tki file)\n"
"    <N>.cfinfo .............. I Constituent info file for right factor\n"
"                                (name taken from tki file)\n"
"    <cond> .................. I Matrix with condensed vectors\n"
"    <uncond> ................ O Matrix with uncondensed vectors\n"
};

static MtxApplication_t *App = NULL;




/* --------------------------------------------------------------------------
   AllocateResult() - Calculate dimension and allocate result matrix

   Description:
     This function cacluates the dimension of M x N, and allocates the
     matrix <UncondMat> which holds the uncondensed vectors.
   -------------------------------------------------------------------------- */

static void AllocateResult()

{
    int i;

    /* Calculate the dimension of M
       ---------------------------- */
    DimM = 0;
    for (i = 0; i < InfoM->nCf; ++i)
	DimM += InfoM->Cf[i].dim * InfoM->Cf[i].mult;

    /* Calculate the dimension of N 
       ---------------------------- */
    DimN = 0;
    for (i = 0; i < InfoN->nCf; ++i)
	DimN += InfoN->Cf[i].dim * InfoN->Cf[i].mult;

    UncondMat = matAlloc(CondMat->field,CondMat->nor,DimM * DimN);
}




/* --------------------------------------------------------------------------
   ReadQMatrices() - Read the Q matrices
   -------------------------------------------------------------------------- */

static void ReadQMatrices()
{
    int i;
    char fn[200];

    for (i = 0; i < TkInfo.nCf; ++i)
    {
	sprintf(fn,"%s.q.%d",TkiName,i+1);
	MTX_LOGD("Reading %s",fn);
	QMat[i] = matLoad(fn);
    }
}


/* --------------------------------------------------------------------------
   Init() - Initialize everything

   Description:
     This function initializes global variables, processes command line
     options and arguments, and reads the .cfinfo and .tki files.

   Arguments:
     <argc>: Number of arguments
     <argv>: Command line arguments

   Return:
     0 = ok, -1 = error
   -------------------------------------------------------------------------- */

static void Init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    appGetArguments(App,3,3);
    TkiName = App->argV[0];

    tkReadInfo(&TkInfo,TkiName);
    InfoM = latLoad(TkInfo.nameM);
    InfoN = latLoad(TkInfo.nameN);
    CondMat = matLoad(App->argV[1]);

    UncondName = App->argV[2];
    AllocateResult();
    ReadQMatrices();
}




/* --------------------------------------------------------------------------
   CalculatePositions() - Calculate basis vector indexes

   Description:
     This function calculates the index of the first basis vector that 
     belongs to the tensor product of two given constituents. The 
     constituents are given by their number <cf> (as defined in the tkinfo
     file), and by a number for both M and N specifying which copy of the
     constituent is meant.

   Arguments:
     <cf>: Constituent index in Tki-File (0-based).
     <num_m>: Used to distiguish the different copies, if the constituent
       appears in M with multiplicity > 1. (0-based)
     <num_n>: Used to distiguish the different copies, if the constituent
       appears in N with multiplicity > 1. (0-based)
     <condpos>: Index of the first basis vector in the condensed module.
     <uncondpos>: Index of the first basis vector in the uncondensed module.
   -------------------------------------------------------------------------- */

static void CalculatePositions(int cf, int num_m, int num_n, int *condpos,
    int *uncondpos)

{
    int cfm = TkInfo.cfIndex[0][cf];
    int cfn = TkInfo.cfIndex[1][cf];
    int startm = 0, startn = 0;
    int i;

    MTX_ASSERT(cfm >= 0 && cfm <= InfoM->nCf);
    MTX_ASSERT(cfn >= 0 && cfn <= InfoN->nCf);
    MTX_ASSERT(InfoM->Cf[cfm].dim == InfoN->Cf[cfn].dim);


    /* Calculate starting position in M x N.
       ------------------------------------- */
    for (i = 0; i < cfm; ++i)
	startm += InfoM->Cf[i].dim * InfoM->Cf[i].mult;
    startm += InfoM->Cf[cfm].dim * num_m;
    for (i = 0; i < cfn; ++i)
	startn += InfoN->Cf[i].dim * InfoN->Cf[i].mult;
    startn += InfoN->Cf[cfn].dim * num_n;
    *uncondpos = startm * DimN + startn;

    /* Calculate starting position in (M x N)e.
       ---------------------------------------- */
    *condpos = 0;
    for (i = 0; i < cf; ++i)
    {
	int m = TkInfo.cfIndex[0][i];
	int n = TkInfo.cfIndex[1][i];
	*condpos += InfoM->Cf[m].mult * InfoN->Cf[n].mult * InfoM->Cf[m].spl;
    }
    *condpos += (num_m * InfoN->Cf[cfn].mult +num_n)* InfoM->Cf[cfm].spl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Uncondense one component of a vector
//  <row>: Number of vector to uncondense (0-based).
//  <cf>: Constituent Index in Tki-File (0-based).

static void UncondenseCf(int row, int cf)
{
   int i;
   int cfm, cfn, mult_m, mult_n, cf_dim;

   MTX_ASSERT(row >= 0);
   MTX_ASSERT(cf >= 0 && cf < TkInfo.nCf);

   cfm = TkInfo.cfIndex[0][cf];     // Constituent index in M
   cfn = TkInfo.cfIndex[1][cf];     // Constituent index in N
   mult_m = InfoM->Cf[cfm].mult;     // Mutiplicity of constituent in M
   mult_n = InfoN->Cf[cfn].mult;     // Mutiplicity of constituent in N
   cf_dim = InfoM->Cf[cfm].dim;      // Dimension of the constituent

   for (i = 0; i < mult_m; ++i) {
      int j;
      for (j = 0; j < mult_n; ++j) {
         int k;
         Matrix_t* condvec;
         int pos, condpos, uncondpos;

         CalculatePositions(cf, i, j, &condpos, &uncondpos);
         condvec = matDupRegion(CondMat, row, condpos, 1, QMat[cf]->nor);
         matMul(condvec, QMat[cf]);
         pos = 0;
         for (k = 1; k <= InfoM->Cf[cfm].dim; ++k) {
            matCopyRegion(UncondMat, row, uncondpos, condvec, 0, pos, 1, cf_dim);
            uncondpos += DimN;
            pos += cf_dim;
         }
         matFree(condvec);
      }
   }
}

/* --------------------------------------------------------------------------
   Uncondense() - Uncondense one vector

   Description:

   Arguments:
     <row>: Number of row to uncondense (starting with 1).
   -------------------------------------------------------------------------- */

static void Uncondense(int row)
{
    int cf;

    /* Loop through all constituents
       ----------------------------- */
    for (cf = 0; cf < TkInfo.nCf; ++cf)
	UncondenseCf(row,cf);
}





/* --------------------------------------------------------------------------
   main() - Program entry point
   -------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
    Init(argc,argv);
    for (int row = 0; row < CondMat->nor; ++row)
	Uncondense(row);
    matSave(UncondMat,UncondName);
    latDestroy(InfoM);
    latDestroy(InfoN);
    if (App != NULL)
	appFree(App);
    return 0;
}



/**
@page prog_tuc tuc - Tensor Uncondense

@see
- @ref prog_precond
- @ref prog_pwkond
- @ref prog_tcond
- @ref prog_zts

@section tuc_syntax Command Line
<pre>
tuc @em Options @em Info @em Cond @em Uncond
</pre>

@par @em Options
Standard options, see @ref prog_stdopts

@par @em Info
Name of tensor condensation data file (.tki file).

@par @em Info
Name of tensor condensation data file (.tki file).

@par @em Cond
Matrix containing condensed vectors.

@par @em Uncond
Output file for uncondensed vectors.


@section tuc_inp Input Files
@par @em Info.tki
Information generated by tcond.
@par @em Cond
Matrix containing condensed vectors.

@section tuc_out Output Files
@par Uncond
Output file for uncondensed vectors.

@section tuc_desc Description
This program is part of the tensor condensation package. It is used to
uncondense one or more vectors, i.e., it calculates the embedding of 
(M⊗N)e<sub>H</sub> into (M⊗N).
The input vectors must be given with respect to the basis defined
by PRECOND, which is also used by TCOND. The uncondensed vectors are 
calculated with respect to the semisimplicity basis produced by PWKOND.

*/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
