////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - This program uncondenses one or more vectors in (M x N)e.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>
#include <string.h>

/* --------------------------------------------------------------------------
   Global data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO

const char *TkiName;		    /* Base name for .tki file */
TkData_t TkInfo;		    /* Data from .tki file */
Lat_Info InfoM, InfoN;		    /* Data from .cfinfo files */
Matrix_t *CondMat;		    /* Condensed vectors (read from file) */
Matrix_t *UncondMat;		    /* Uncondensed vectors (calculated here) */
const char *UncondName;		    /* Name of output file */
long DimM, DimN;		    /* Dimensions of M and N */
Matrix_t *QMat[LAT_MAXCF];	    /* Q matrices */

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
    for (i = 0; i < InfoM.NCf; ++i)
	DimM += InfoM.Cf[i].dim * InfoM.Cf[i].mult;

    /* Calculate the dimension of N 
       ---------------------------- */
    DimN = 0;
    for (i = 0; i < InfoN.NCf; ++i)
	DimN += InfoN.Cf[i].dim * InfoN.Cf[i].mult;

    UncondMat = MatAlloc(CondMat->Field,CondMat->Nor,DimM * DimN);
}




/* --------------------------------------------------------------------------
   ReadQMatrices() - Read the Q matrices
   -------------------------------------------------------------------------- */

static void ReadQMatrices()
{
    int i;
    char fn[200];

    for (i = 0; i < TkInfo.NCf; ++i)
    {
	sprintf(fn,"%s.q.%d",TkiName,i+1);
	MESSAGE(1,("Reading %s\n",fn));
	QMat[i] = MatLoad(fn);
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

static int Init(int argc, const char **argv)

{
    /* Initialize the MeatAxe library, process command line
       ---------------------------------------------------- */
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (AppGetArguments(App,3,3) != 3)
	return -1;

    /* Read input files
       ---------------- */
    TkiName = App->ArgV[0];
    if (TK_ReadInfo(&TkInfo,TkiName) != 0)
	return -1;
    if (Lat_ReadInfo(&InfoM,TkInfo.NameM) != 0 ||
        Lat_ReadInfo(&InfoN,TkInfo.NameN) != 0)
	return -1;
    CondMat = MatLoad(App->ArgV[1]);
    if (CondMat == NULL)
	return -1;

    /* Other initializations
       --------------------- */
    UncondName = App->ArgV[2];
    AllocateResult();
    ReadQMatrices();
    return 0;
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
    int cfm = TkInfo.CfIndex[0][cf];
    int cfn = TkInfo.CfIndex[1][cf];
    int startm = 0, startn = 0;
    int i;

    MTX_VERIFY(cfm >= 0 && cfm <= InfoM.NCf);
    MTX_VERIFY(cfn >= 0 && cfn <= InfoN.NCf);
    MTX_VERIFY(InfoM.Cf[cfm].dim == InfoN.Cf[cfn].dim);


    /* Calculate starting position in M x N.
       ------------------------------------- */
    for (i = 0; i < cfm; ++i)
	startm += InfoM.Cf[i].dim * InfoM.Cf[i].mult;
    startm += InfoM.Cf[cfm].dim * num_m;
    for (i = 0; i < cfn; ++i)
	startn += InfoN.Cf[i].dim * InfoN.Cf[i].mult;
    startn += InfoN.Cf[cfn].dim * num_n;
    *uncondpos = startm * DimN + startn;

    /* Calculate starting position in (M x N)e.
       ---------------------------------------- */
    *condpos = 0;
    for (i = 0; i < cf; ++i)
    {
	int m = TkInfo.CfIndex[0][i];
	int n = TkInfo.CfIndex[1][i];
	*condpos += InfoM.Cf[m].mult * InfoN.Cf[n].mult * InfoM.Cf[m].spl;
    }
    *condpos += (num_m * InfoN.Cf[cfn].mult +num_n)* InfoM.Cf[cfm].spl;
}


/* --------------------------------------------------------------------------
   UncondenseCf() - Uncondense one component of a vector

   Description:

   Arguments:
     <row>: Number of vector to uncondense (0-based).
     <cf>: Constituent Index in Tki-File (0-based).
   -------------------------------------------------------------------------- */

static void UncondenseCf(int row, int cf)
{
    int i;
    int cfm, cfn, mult_m, mult_n, cf_dim;

    MTX_VERIFY(row >= 0);
    MTX_VERIFY(cf >= 0 && cf < TkInfo.NCf);

    cfm = TkInfo.CfIndex[0][cf];    /* Constituent index in M */
    cfn = TkInfo.CfIndex[1][cf];    /* Constituent index in N */
    mult_m = InfoM.Cf[cfm].mult;    /* Mutiplicity of constituent in M */
    mult_n = InfoN.Cf[cfn].mult;    /* Mutiplicity of constituent in N */
    cf_dim = InfoM.Cf[cfm].dim;	    /* Dimension of the constituent */

    for (i = 0; i < mult_m; ++i)
    {
	int j;
	for (j = 0; j < mult_n; ++j)
	{
	    int k;
	    Matrix_t *condvec;
	    int pos, condpos, uncondpos;

	    CalculatePositions(cf,i,j,&condpos,&uncondpos);
	    condvec = MatCut(CondMat,row,condpos,1,QMat[cf]->Nor);
	    MatMul(condvec,QMat[cf]);
	    pos = 0;
	    for (k = 1; k <= InfoM.Cf[cfm].dim; ++k)
	    {
		MatCopyRegion(UncondMat,row,uncondpos,condvec,0,pos,1,
		    cf_dim);
		uncondpos += DimN;
		pos += cf_dim;
	    }
	    MatFree(condvec);
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
    for (cf = 0; cf < TkInfo.NCf; ++cf)
	UncondenseCf(row,cf);
}





/* --------------------------------------------------------------------------
   main() - Program entry point
   -------------------------------------------------------------------------- */

int main(int argc, const char **argv)

{
    int row;

    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return 1;
    }
    for (row = 0; row < CondMat->Nor; ++row)
	Uncondense(row);
    MatSave(UncondMat,UncondName);
    if (App != NULL)
	AppFree(App);
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

