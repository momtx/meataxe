////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Calculate the incidences between mountains
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>
#include <stdlib.h>

MatRep_t *Rep;				// Generators
Matrix_t *bild[LAT_MAXCF];		// Image of peak word (gkond)
int nmount = 0;				// Number of mountains
Matrix_t * mountlist[MAXCYCL];		// List of mountains
int MountDim[MAXCYCL];			// Dim. of mountains
Matrix_t **proj[MAXCYCL];		// Projections of mountains
int moffset[LAT_MAXCF];			// Number of first mountain
int *Class[MAXCYCL];			// Classes of vectors
BitString_t *subof[MAXCYCL];		// Incidence matrix
static Lat_Info LI;			// Data from .cfinfo file

static MtxApplicationInfo_t AppInfo = { 
"mkinc", "Mountains and incidence matrix",
"\n"
"SYNTAX\n"
"    mkinc [<Options>] <Name>\n"
"\n"
"ARGUMENTS\n"
"    <Name> .................. Name of the representation\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -G ...................... GAP output (implies -Q)\n"
"\n"
"FILES\n"
"    <Name>.cfinfo ........... IO Constituent info file\n"
"    <Name><Cf>.{1,2...} ..... I Generators on the constituents\n"
"    <Name><Cf>.{1,2...}k .... I Generators on the condensed modules\n"
"    <Name><Cf>.v ............ I Cyclic submodules, generated by MKCYCL\n"
"    <Name><Cf>.im ........... I Images used for condensation\n"
"    <Name><Cf>.k ............ I Uncondense matrices\n"
"    <Name>.v ................ O Mountains\n"
"    <Name>.mnt .............. O Mountain dimensions and classes of cyclic\n"
"                                submodules corresponding to the mountains\n"
"    <Name>.inc .............. O Incidence matrix\n"
};

static MtxApplication_t *App = NULL;

int opt_G = 0;	    // GAP output


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Read generators and images of peak words

static void readFiles(const char *basename)
{
    latReadInfo(&LI,basename);
    Rep = mrLoad(LI.BaseName,LI.NGen);

    // Load the .im matrices
    for (int i = 0; i < LI.nCf; ++i)
    {
        char fn[200];
	sprintf(fn,"%s%s.im",LI.BaseName,latCfName(&LI,i));
	bild[i] = matLoad(fn);
	matEchelonize(bild[i]);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    opt_G = appGetOption(App,"-G --gap");
//    if (opt_G) 
//	MtxMessageLevel = -100;
    appGetArguments(App,1,1);
    MTX_LOGI("Start mkinc - Find mountains and their incidence relation");
    readFiles(App->argV[0]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Write mountains, dimensions and classes (files XXX.v and XXX.mnt)

static void writeMountains()
{
    char fn[200];
    int i;

    // Write dimensions and classes
    FILE* f = sysFopen(strcat(strcpy(fn,LI.BaseName),".mnt"),"w");
    MTX_LOGD("Writing dimensions and classes to %s",fn);
    for (i = 0; i < nmount; ++i)
    {
	int *p;
	fprintf(f,"%4d %4d  %d ",i,MountDim[i],Class[i][0]);
	for (p = Class[i] + 1; *p >= 0; ++p)
	    fprintf(f,"%d ",*p + 1);
	fprintf(f,"-1\n");
    }
    fclose(f);

    // Write mountains
    strcat(strcpy(fn,LI.BaseName),".v");
    MTX_LOGD("Writing mountains to %s",fn);
    MtxFile_t* mountainsFile = mfCreate(fn,ffOrder,nmount,Rep->Gen[0]->noc);
    for (i = 0; i < nmount; ++i)
    {
	ffWriteRows(mountainsFile, mountlist[i]->data, 1, mountlist[i]->noc);
	matFree(mountlist[i]);	// We don't need them for step 2*
        mountlist[i] = NULL;
    }
    mfClose(mountainsFile);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Takes a vector, spins it up and checks if it is a new mountain. If yes, it is added to the
// mountain list and projection on each condensed module is calculated.
// Returns 1 if a new mountain has been found, 0 if not.

static int newMountain(Matrix_t *vec, int cf)

{
    Matrix_t *span, *backproj;
    int i;

    // Spin up the vector and project back onto the condensed module where it came from.
    span = SpinUp(vec,Rep,SF_FIRST|SF_SUB,NULL,NULL);
    MTX_LOG2("Next vector spins up to %d",span->nor);
    backproj = QProjection(bild[cf],span);
    matEchelonize(backproj);

    // Check if it is a new mountain
    for (i = moffset[cf]; i < nmount; ++i)
    {
	if (backproj->nor == proj[i][cf]->nor)
	{
	    if (IsSubspace(proj[i][cf],backproj,0))
		break;
	}
    }

    // If it is new, add it to the list and calculate the other projections.
    // Otherwise just forget it.
    if (i >= nmount)
    {
	int k;

	if (nmount >= MAXCYCL)
	    mtxAbort(MTX_HERE,"TOO MANY MOUNTAINS, INCREASE MAXCYCL");
	proj[nmount] = NALLOC(Matrix_t *,LI.nCf);

    	MTX_LOG2("New Mountain %d",(int)i);
	for (k = 0; k < LI.nCf; ++k)
	{
    	    MTX_LOG2("Projecting on %d",k);
	    if (k == cf)
	    	proj[nmount][cf] = backproj;
	    else
	    {
		proj[nmount][k] = QProjection(bild[k],span);
		matEchelonize(proj[nmount][k]);
	    }
	}
	mountlist[nmount] = vec;
	MountDim[nmount] = span->nor;
	++nmount;
	matFree(span);
	return 1;
    }
    else
    {
	matFree(backproj);
	matFree(span);
	matFree(vec);
	return 0;
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////// 

// Finds all cyclic subspaces of the projection of a mountain onto `its' irreducible.

static void makeclass(int mnt, int cf, Matrix_t* vectors)
{
   char* tmp = NALLOC(char, vectors->nor + 2);
   size_t nvec = 0;
   MTX_LOG2("Making equivalence class");
   for (int k = 0; k < vectors->nor; ++k) {
      Matrix_t* vec = matCutRows(vectors, k, 1);
      tmp[k] = 0;
      if (IsSubspace(vec, proj[mnt][cf], 1)) {
         tmp[k] = 1;
         ++nvec;
      }
      matFree(vec);
   }

   int* p = Class[mnt] = NALLOC(int, nvec + 2);
   *p++ = nvec;
   for (int k = 0; k < vectors->nor; ++k) {
      if (tmp[k]) {
         *p++ = k;
      }
   }
   *p = -1;
   sysFree(tmp);
}


//////////////////////////////////////////////////////////////////////////////////////////////////// 

// Make all mountains and calculate the projections of mountains on condensed modules.

static void findMountains()
{
   Matrix_t* vectors, * vec, * uk;
   char fn[200];
   int cf;
   int i;

   MTX_LOGI("Step 1 (Mountains)");
   nmount = 0;
   for (cf = 0; cf < LI.nCf; ++cf) {
      // Read the vectors and the uncondense matrix
      sprintf(fn, "%s%s.v", LI.BaseName, latCfName(&LI, cf));
      vectors = matLoad(fn);
      sprintf(fn, "%s%s.k", LI.BaseName, latCfName(&LI, cf));
      uk = matLoad(fn);

      // Try each vector
      moffset[cf] = nmount;
      for (i = 0; i < vectors->nor; ++i) {
         vec = matCutRows(vectors, i, 1);
         matMul(vec, uk);       // Uncondense
         if (newMountain(vec, cf)) {
            makeclass(nmount - 1, cf, vectors);
         }
      }
      LI.Cf[cf].nmount = nmount - moffset[cf];

      matFree(vectors);
      matFree(uk);
      MTX_LOGI("%s%s: %ld mountain%s", LI.BaseName, latCfName(&LI, cf),
            LI.Cf[cf].nmount, LI.Cf[cf].nmount != 1 ? "s" : "");

   }
   MTX_LOGI("Total: %d mountain%s", nmount, nmount != 1 ? "s" : "");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void writeIncidenceMatrix()
{
   MtxFile_t* file = mfOpen(strEprintf("%s.inc", LI.BaseName), "wb");
   MTX_LOGD("Writing incidence matrix (%s)", file->name);
   uint32_t l = nmount;
   mfWrite32(file, &l, 1);
   for (int i = 0; i < nmount; ++i) {
      bsWrite(subof[i], file);
   }
   mfClose(file);

   latWriteInfo(&LI);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void calculateIncidences()
{
    int i, k;
    int cfi, cfk;	// Comp. factor corresponding to mountain i,j

    MTX_LOGI("Step 2 (Incidences)");

    // Allocate memory for the incidence matrix
    for (i = 0; i < nmount; ++i)
	subof[i] = bsAlloc(nmount);

    // Calculate the incidences
    for (cfi = i = 0; i < nmount; ++i)
    {
	if (i == moffset[cfi])
	    MTX_LOGI("%s%s",LI.BaseName,latCfName(&LI,cfi));
	for (cfk=0, k = 0; k < nmount; ++k)
	{
	    int isubk, ksubi;
	    ksubi = IsSubspace(proj[k][cfk],proj[i][cfk],0);
	    isubk = IsSubspace(proj[i][cfi],proj[k][cfi],0);
	    if (ksubi < 0 || isubk < 0)
	    {
		mtxAbort(MTX_HERE,"Subspace comparison failed");
		return;
	    }
	    if (ksubi)
		bsSet(subof[k],i);
	    if (isubk)
		bsSet(subof[i],k);
	    if (cfk < LI.nCf && k == moffset[cfk+1]-1)
		++cfk;
    	}
    	if (cfi < LI.nCf && i == moffset[cfi+1]-1)
	    ++cfi;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
   
/// Writes the mountain list and the incidence matrix in GAP format to stdout.

static void WriteResultGAP()
{
    int i,  j;
    int *p;

    // Write the incidence matrix
    printf("MeatAxe.NMount := %d;\n", nmount);
    printf("MeatAxe.Incidences := [\n");
    for (i = 0; i < nmount; ++i)
    {
	printf("BlistList([" );
	for (j = 0 ; j < nmount ; j++)
	{
      	    printf(j < (nmount - 1) ? "%s," : "%s], [1])", 
		bsTest(subof[i],j) ? "1" : "0" ) ;
	}

	printf( i < nmount-1 ? ",\n" : "] ;\n" ) ;
    }

    // Write dimensions and classes
    printf("MeatAxe.Dimensions := [");
    for (i = 0; i < nmount-1 ; ++i)
	printf("%d,", MountDim[i]);
    printf("%d] ;\n", MountDim[nmount-1]);

    printf( "MeatAxe.Classes := [\n");
    for (i = 0; i < nmount; ++i)
    {
        printf("[%d",Class[i][0]);
	for (p = Class[i] + 1; *p >= 0 ; ++p)
	{
	    printf(",%d", *p + 1);
	}
	printf("]");
	printf(i < nmount-1 ? ",\n" : "] ;\n");
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void cleanup()
{
   for (int m = 0; m < nmount; ++m) {
      for (int k = 0; k < LI.nCf; ++k) {
         matFree(proj[m][k]);
      }
   }
   for (int i = 0; i < LI.nCf; ++i) {
      matFree(bild[i]);
   }
   for (int i = 0; i < nmount; ++i) {
      bsFree(subof[i]);
   }
   mrFree(Rep);
   latCleanup(&LI);
   appFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    init(argc,argv);
    findMountains();
    writeMountains();
    calculateIncidences();
    writeIncidenceMatrix();
    if (opt_G) 
	WriteResultGAP();
    cleanup();
    return 0;
}


/**
@page prog_mkinc mkinc - Find Mountains

@section mkinc_syntax Command Line
<pre>
mkinc [@em Options] [-G] @em Name
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -G
  Produce output in GAP format.
@par @em Name
  Name of the representation.

@section mkinc_inp Input Files
@par @em Name.cfinfo
  Constituent info file.
@par @em NameCf.1, @em NameCf.2, ...
  Generators on the irreducible constituents.
@par @em NameCf.1k, @em NameCf.2k, ...
  Generators on the corresponding condensed modules.
@par @em NameCf.v
  Cyclic submodules, generated by @ref prog_mkcycl "mkcycl".
@par @em NameCf.im
  Peak word images.
@par @em NameCf.k
  Uncondense matrices.

@section mkinc_out Output Files
@par @em Name.v
  Mountains.
@par @em Name.mnt
  Mountain dimensions and classes of cyclic submodules corresponding to the mountains.
@par @em Name.inc
  Incidence matrix.

@section mkinc_desc Description
The MKINC program is part of the @ref sec_progs_lattice "Submodule Lattice Package".
It is invoked after @ref prog_pwkond "pwkond" has calculated the cyclic submodules.
MKCYCL runs in two steps. During the first step, all cyclic submodules are uncondensed,
giving the local submodules, the "mountains", of the original module. Then, each
local submodule is projected back to the condensed module, and all
cyclic vectors which are contained in the image are found.
At the end of this step, there is a list of local submodules and,
for each local submodule, a list of cyclic subspaces in the
condensed module. This information is written to @em Name.mnt.

In the second step, @b mkinc computes the incidence relation between
local submodules. The result is a matrix which contains a 1 for
each incidence. This matrix is written to the file @em Name.inc.


@section mkinc_impl Implementation Details
The whole calculation of step 2 is done in the condensed modules. 
This is possible because incidences between local submodules do not 
change if they are condensed. Usually this saves a lot of both memory 
and CPU time because one does not have to keep all mountains 
simultaneously, and the condensed modules have a smaller dimension.
**/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
