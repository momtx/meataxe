////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Calculates the socle series of a module.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"



/* --------------------------------------------------------------------------
   Global data
   -------------------------------------------------------------------------- */


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



////////////////////////////////////////////////////////////////////////////////////////////////////

static void ReadConstituents()
{
    int i;
    char fn[200];

    for (i = 0; i < LI.nCf; ++i)
    {
	/* Read generators
	   --------------- */
	sprintf(fn,"%s%s.std", LI.BaseName,latCfName(&LI,i));
	CfRep[i] = mrLoad(fn,LI.NGen);

	/* Read spinup script for standard basis
	   ------------------------------------- */
	sprintf(fn,"%s%s.op",LI.BaseName,latCfName(&LI,i));
	OpTableMat[i] = imatLoad(fn);
	if (ConvertSpinUpScript(OpTableMat[i]))
	{
	    static int first_time = 1;
	    if (first_time)
	    {
		MTX_LOGD("Converting spinup script from MeatAxe 2.3 format");
		first_time = 0;
	    }
	}

	/* Read peak word kernel.
	   ---------------------- */
	sprintf(fn,"%s%s.k",LI.BaseName,latCfName(&LI,i));
	MTX_LOGD("Taking seed vectors from %s",fn);
	seed[i] = matLoad(fn);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    MaxLen = appGetIntOption(App,"-l --max-length",0,0,1000);
    appGetArguments(App,1,1);

    latReadInfo(&LI,App->argV[0]);
    Rep = mrLoad(App->argV[0],LI.NGen);
    ReadConstituents();
    Dimension = Rep->Gen[0]->nor;

    WGen = wgAlloc(Rep);
    LI.NSocles = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void WriteBasis(const Matrix_t *basis)
{
    char name[200];
    sprintf(name, "%s.soc",App->argV[0]);
    MTX_LOGD("Writing basis to %s",name);
    matSave(basis,name);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int NextLayer()
{
    int flag = 0;

    SocDim = 0;
    Matrix_t* bas = matAlloc(ffOrder,Dimension,Dimension);
    int* cfvec = NALLOC(int,LI.nCf);

    for (int j = 0; j < LI.nCf; j++)
    {
        if (LI.Cf[j].peakWord == 0)
            mtxAbort(MTX_HERE,"Missing peak word for constituent %d - run pwkond!",j);
	
	// determine a basis for the corresponding part in the socle
	Matrix_t *partbas;
	if (seed[j]->nor != 0) 
	{
	    const int dimEndoS = LI.Cf[j].spl;
	    partbas = HomogeneousPart(Rep,CfRep[j],seed[j],OpTableMat[j],dimEndoS);
	    if (partbas == NULL)
		mtxAbort(MTX_HERE,"The module is too big");
	}
	else 
	    partbas = matDup(seed[j]);
                
	// Append the new basis to the old one.
	cfvec[j] = partbas->nor / LI.Cf[j].dim;
	matCopyRegion(bas,SocDim,0,partbas,0,0,-1,-1);
	SocDim += partbas->nor;
	MTX_LOG2("Socle dimension of %s is %d", latCfName(&LI,j), partbas->nor);
        matFree(partbas);
    }


    // makes the output
    ++SocLen;
    MTX_XLOGI(msg) {
       sbPrintf(msg, "Socle %d: %d =",SocLen,SocDim);
       flag = 0;
       for (int j = 0; j < LI.nCf; j++)
       {
          if (cfvec[j] <= 0) 
             continue;
          if (flag++ > 0)
             sbPrintf(msg," +");
          if (cfvec[j] == 1)
             sbPrintf(msg, " %s",latCfName(&LI,j));
          else
             sbPrintf(msg, " %d*%s",cfvec[j],latCfName(&LI,j));
       }
    }
    latAddSocle(&LI,cfvec);

  
    sysFree(cfvec);


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
	    mat = matCutRows(basis,basis->nor - SocDim,SocDim);
	    matMul(bas,mat);
	    matFree(mat);
	    matCopyRegion(basis,basis->nor - SocDim,0,bas,0,0,SocDim,-1);
	    WriteBasis(basis);
	}
        matFree(bas);
        return 1;
    }

    // Extend the basis of the socle to a basis of the whole module.
    Matrix_t* echbas = matAlloc(bas->field,bas->noc,bas->noc);
    matCopyRegion(echbas,0,0,bas,0,0,-1,-1);
    matEchelonize(bas);
    for (uint32_t i = bas->nor; i < bas->noc; ++i)
	ffInsert(matGetPtr(echbas,i),bas->pivotTable[i],FF_ONE);
    matFree(bas);
    bas = echbas;

    // multiplying the last two basis changes
    if (basis == NULL)
	basis = matDup(bas);
    else
    {
	Matrix_t* mat = matCutRows(basis,basis->nor - Dimension,Dimension);
	Matrix_t* stgen = matDup(bas);
	matMul(stgen, mat);
	matCopyRegion(basis,basis->nor - Dimension,0,stgen,0,0,Dimension,-1);
	matFree(stgen);
	matFree(mat);
    }

    // exiting when the first MaxLen socles have been calculated
    if (SocLen == MaxLen)
    {
	WriteBasis(basis);
        matFree(bas);
	return 1;
    }

    // factorizing with the socle

    Matrix_t* basi = matInverse(bas);

    // the kernels
    for (int j = 0; j < LI.nCf; j++)
    {
	Matrix_t *tmp;
	if (seed[j] == NULL)
	    continue;
	matMul(seed[j], basi);
	tmp = matCut(seed[j],0,SocDim,-1,-1);
	matFree(seed[j]);
	seed[j] = tmp;
	matEchelonize(seed[j]);
    }

    // the generators
    for (int i = 0; i < LI.NGen; i++)
    {
        Matrix_t *stgen = matDup(bas);
        matMul(stgen,Rep->Gen[i]);
        matMul(stgen,basi);
        matFree(Rep->Gen[i]);
	Rep->Gen[i] = matCut(stgen,SocDim,SocDim,-1,-1);
        matFree(stgen);
    }
    Dimension = Rep->Gen[0]->nor;
    MTX_LOGD("Reduced to dimension %d",Dimension);

    matFree(bas);
    matFree(basi);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void cleanup()
{
   if (basis != NULL)
      matFree(basis);
   wgFree(WGen);
   for (int i = 0; i < LI.nCf; ++i) {
      mrFree(CfRep[i]);
      imatFree(OpTableMat[i]);
      matFree(seed[i]);
   }
   latCleanup(&LI);
   mrFree(Rep);
   appFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char **argv)
{
    init(argc,argv);

    while (NextLayer() == 0);

    latWriteInfo(&LI);
    if (SocDim != Dimension)
    {
	MTX_LOGI("Warning: Calculation aborted at dimension %d of %d", SocDim,Dimension);
    }
    cleanup();
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
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
