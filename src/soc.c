////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Calculates the socle series of a module.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"



/* --------------------------------------------------------------------------
   Global data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO

static MtxApplicationInfo_t AppInfo = { 
"soc", "Socle series", 
"\n"
"SYNTAX\n"
"    soc " MTX_COMMON_OPTIONS_SYNTAX " [-l <length>] <Name>\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -l <length> .............. Calculate only the first <length> layers\n"
"\n"
"FILES\n"
"    <Name>.cfinfo ........... IO Constituent info file\n"
"    <Name><Cf>.std.{1,2...} . I  Generators on constituents\n"
"    <Name><Cf>.op ........... I  Spin-up script for standard basis\n"
"    <Name><Cf>.k ............ I  Uncondense matrix\n"
"    <Name>.soc .............. O  Matrix for basis change\n"
};

static MtxApplication_t *App = NULL;

static int MaxLen = 0;			/* Max. length of s.series (0 = no limit) */
static int SocLen = 0;			/* Socle length */
static Lat_Info LI;			/* Data from .cfinfo file */
static IntMatrix_t *OpTableMat[LAT_MAXCF];  /* Operations, written to <Op> */
static MatRep_t *Rep;			/* The module */
static int Dimension = 0;		/* Module dimension */
static MatRep_t *CfRep[LAT_MAXCF];	/* Constituents in standard basis */
static Matrix_t *seed[LAT_MAXCF];	/* Kernel of the peak words */
static WgData_t *WGen;
static int SocDim;
static Matrix_t *basis = NULL;		/* Basis corresponding to Loewy series */



/* --------------------------------------------------------------------------
   ParseCommandLine() - Process command line options and arguments
   -------------------------------------------------------------------------- */

static int ParseCommandLine()

{
    MaxLen = AppGetIntOption(App,"-l --max-length",0,0,1000);
    if (AppGetArguments(App,1,1) != 1)
	return -1;
    return 0;
}


/* --------------------------------------------------------------------------
   ReadConstituents() - Read generators and optables for all constituents
   -------------------------------------------------------------------------- */

static int ReadConstituents()

{

    int i;
    char fn[200];

    for (i = 0; i < LI.NCf; ++i)
    {
	/* Read generators
	   --------------- */
	sprintf(fn,"%s%s.std", LI.BaseName,Lat_CfName(&LI,i));
	CfRep[i] = MrLoad(fn,LI.NGen);
	if (CfRep[i] == NULL)
	    return -1;

	/* Read spinup script for standard basis
	   ------------------------------------- */
	sprintf(fn,"%s%s.op",LI.BaseName,Lat_CfName(&LI,i));
	OpTableMat[i] = ImatLoad(fn);
	if (ConvertSpinUpScript(OpTableMat[i]))
	{
	    static int first_time = 1;
	    if (first_time)
	    {
		MESSAGE(1,("Converting spinup script from MeatAxe 2.3 format\n"));
		first_time = 0;
	    }
	}

	/* Read peak word kernel.
	   ---------------------- */
	sprintf(fn,"%s%s.k",LI.BaseName,Lat_CfName(&LI,i));
	MESSAGE(1,("Taking seed vectors from %s\n",fn));
	seed[i] = MatLoad(fn);
    }
    return 0;
}




/* --------------------------------------------------------------------------
   Init() - Program initialization
   -------------------------------------------------------------------------- */

static int Init(int argc, const char **argv)

{
    /* Process the command line
       ------------------------ */
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (ParseCommandLine() != 0)
    {
	MTX_ERROR("Error in command line");
	return -1;
    }

    /* Read input files
       ---------------- */
    if (Lat_ReadInfo(&LI,App->ArgV[0]) != 0)
	return -1;
    if ((Rep = MrLoad(App->ArgV[0],LI.NGen)) == NULL)
	return -1;
    if (ReadConstituents() != 0)
	return -1;
    Dimension = Rep->Gen[0]->Nor;

    WGen = WgAlloc(Rep);
    LI.NSocles = 0;

    return 0;
}


static void WriteBasis(const Matrix_t *basis)

{
    char name[200];
    sprintf(name, "%s.soc",App->ArgV[0]);
    MESSAGE(1,("Writing basis to %s\n",name));
    MatSave(basis,name);
}



static int NextLayer()

{
    int *cfvec;
    Matrix_t *bas, *echbas;
    int i, j;
    int flag;
    Matrix_t *basi;

    SocDim = 0;
    if ((bas = MatAlloc(FfOrder,Dimension,Dimension)) == NULL)
    {
	MTX_ERROR("Cannot allocate basis matrix");
	return -1;
    }

    if ((cfvec = NALLOC(int,LI.NCf)) == NULL)
    {
	MTX_ERROR("Cannot allocate vector table");
	return -1;
    }


    for (j = 0; j < LI.NCf; j++)
    {
	Matrix_t *sed;
	int dimendS;
	Matrix_t *partbas;

        if (LI.Cf[j].peakword <= 0)
	{
            MTX_ERROR1("Missing peak word for constituent %d - run pwkond!",j);
	    return -1;
	}
	
	sed = seed[j];
	dimendS = LI.Cf[j].spl;


	/* determines a basis for the corresponding part in the socle
	   ---------------------------------------------------------- */
	if (sed->Nor != 0) 
	{
	    partbas = HomogeneousPart(Rep,CfRep[j],sed,OpTableMat[j],dimendS);
	    if (partbas == NULL)
	    {
		MTX_ERROR("The module is too big");
		return -1;
	    }
	}
	else 
	    partbas = MatDup(sed);
                
	/* Append the new basis to the old one.
	   ------------------------------------ */
	cfvec[j] = partbas->Nor / LI.Cf[j].dim;
	MatCopyRegion(bas,SocDim,0,partbas,0,0,-1,-1);
	SocDim += partbas->Nor;
	MESSAGE(2,("Socle dimension of %s is %d\n", Lat_CfName(&LI,j),
	    partbas->Nor));
    }


/* ----------------
   makes the output
   ---------------- */
    ++SocLen;
    MESSAGE(0,("Socle %d: %d =",SocLen,SocDim));
    flag = 0;
    for (j = 0; j < LI.NCf; j++)
    {
	if (cfvec[j] <= 0) 
	    continue;
	if (flag++ > 0)
	    MESSAGE(0,(" +"));
	if (cfvec[j] == 1)
	    MESSAGE(0,(" %s",Lat_CfName(&LI,j)));
	else
	    MESSAGE(0,(" %d*%s",cfvec[j],Lat_CfName(&LI,j)));
    }
    MESSAGE(0,("\n"));
    Lat_AddSocle(&LI,cfvec);

  
    SysFree(cfvec);


/* --------------------------------------------------
   exiting when the socle is already the whole module
   -------------------------------------------------- */

    if (SocDim == Dimension)
    { 
	/* multiplies the last two basistransformations
	   -------------------------------------------- */
	if (basis == NULL)
	{
	    WriteBasis(bas);
	}
	else
	{
	    Matrix_t *mat;
	    mat = MatCutRows(basis,basis->Nor - SocDim,SocDim);
	    MatMul(bas,mat);
	    MatFree(mat);
	    MatCopyRegion(basis,basis->Nor - SocDim,0,bas,0,0,SocDim,-1);
	    WriteBasis(basis);
	}
        return 1;
    }

    
    
    /* Extend the basis of the socle to a basis of the whole module.
       ------------------------------------------------------------- */
    echbas = MatAlloc(bas->Field,bas->Noc,bas->Noc);
    MatCopyRegion(echbas,0,0,bas,0,0,-1,-1);
    MatEchelonize(bas);
    for (i = bas->Nor; i < bas->Noc; ++i)
	FfInsert(MatGetPtr(echbas,i),bas->PivotTable[i],FF_ONE);
    MatFree(bas);
    bas = echbas;

	    

/* multiplying the last two basischanges
   ------------------------------------- */
    if (basis == NULL)
	basis = MatDup(bas);
    else
    {
	Matrix_t *mat, *stgen;

	mat = MatCutRows(basis,basis->Nor - Dimension,Dimension);
	stgen = MatDup(bas);
	MatMul(stgen, mat);
	MatCopyRegion(basis,basis->Nor - Dimension,0,stgen,0,0,Dimension,-1);
	MatFree(mat);
	MatFree(stgen);
    }

/* ------------------------------------------------------
   exiting when the first MaxLen socles have been calculated
   ------------------------------------------------------ */

    if (SocLen == MaxLen)
    {
	WriteBasis(basis);
	return 1;
    }

/* --------------------------
   factorizing with the socle
   -------------------------- */

    basi = MatInverse(bas);

/* the kernels
   ----------- */
    for (j = 0; j < LI.NCf; j++)
    {
	Matrix_t *tmp;
	if (seed[j] == NULL)
	    continue;
	MatMul(seed[j], basi);
	tmp = MatCut(seed[j],0,SocDim,-1,-1);
	MatFree(seed[j]);
	seed[j] = tmp;
	MatEchelonize(seed[j]);
    }

    /* the generators
       -------------- */
    for (i = 0; i < LI.NGen; i++)
    {
        Matrix_t *stgen = MatDup(bas);
        MatMul(stgen,Rep->Gen[i]);
        MatMul(stgen,basi);
        MatFree(Rep->Gen[i]);
	Rep->Gen[i] = MatCut(stgen,SocDim,SocDim,-1,-1);
        MatFree(stgen);
    }
    Dimension = Rep->Gen[0]->Nor;
    MESSAGE(1,("Reduced to dimension %d\n",Dimension));

    MatFree(bas);
    MatFree(basi);

    return 0;
}




int main( int argc, const char **argv)
{
    int rc;

    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Program initialization failed");
	return 1;
    }


    while ((rc = NextLayer()) == 0);
    if (rc < 0)
    {
	MTX_ERROR("Error calculating socle series");
	return 1;
    }

    /* Clean up
       -------- */
    Lat_WriteInfo(&LI);
    WgFree(WGen);
    if (SocDim != Dimension)
    {
	MESSAGE(0,("Warning: Calculation aborted at dimension %d of %d\n",
	    SocDim,Dimension));
    }
    if (App != NULL) AppFree(App);
    return 0;
}

    

/**
@page prog_soc soc - Socle Series

@section soc_syntax Command Line
<pre>
soc @em Options [-l @em MaxLength] @em Module
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -l @em MaxLength
  Maximum number of layers to compute.
@par @em Module
  Module name.

@section soc_inp Input Files
@par @em Name.cfinfo.
  Constituent information.
@par @em NameCf.std.1, @em NameCf.std.2, ...
  Generators on the irreducible constituents.
@par @em NameCf.op
  Spin-up script for the standard basis.
@par @em NameCf.k
  Uncondense matrix

@section soc_out Output Files
@par @em Name.cfinfo.
  Socle information, see description.
@par @em Name.soc.
  A basis reflecting the Loewy structure.

@section soc_desc Description
This program determines the Loewy structure of a module by calculating the socles.
Before using the program, you must run @ref prog_chop "chop" and @ref prog_pwkond "pwkond"
with the "-t" option. For example,
<pre>
chop m11
pwkond -t m11
soc m11
</pre>
For each layer of the socle series, the program prints the dimension
and the multiplicities of the irreducible constituents in this layer.
This information is also written to the cfinfo file. The following 
example shows the relevant portion of the cfinfo file:
<pre>
CFInfo.NSocles := 5;
CFInfo.Socles := [[1,0,0],[0,1,1],[2,0,0],[0,1,1],[1,0,0]];
</pre>
The numbers in @c CFInfo.socles are the multiplicities of the
irreducible constituents for each layer of the socle series.

Using the "-l" option, you can specify a maximum length.  After @em MaxLength socles
have been calculated, the program prints a warning and stops. 

A basis basis reflecting the Loewy structure of the module is 
written to @em Name.soc.
Note: @b soc always writes a basis of the full space.
If the socle series is not calculated completely because the maximum length
has been reached, the partial basis found so far is extended with random
vectors to form a complete basis.

@section soc_impl Implementation Details
This program uses an algorithm by Magdolna Szöke, see @ref Sz98 "[Sz98]".

**/
