////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Uncondense vectors
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <stdlib.h>
#include <string.h>


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */



static MtxApplicationInfo_t AppInfo = { 
"zuk", "Uncondense Vectors", 
"SYNTAX\n"
"    zuk " MTX_COMMON_OPTIONS_SYNTAX " <Vectors> <Orbits> <Result>\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"\n"
"FILES\n"
"    <Vectors> ............... I A matrix (#columns = #orbits)\n"
"    <Orbits> ................ I Orbits file, produces by ZMO\n"
"    <Result> ................ O Uncondensed vectors\n"
};

static MtxApplication_t *app = NULL;
static const char *fileNameInp;
static const char *fileNameOrbits;
static const char *fileNameOut;
static IntMatrix_t *orbits = NULL;
static IntMatrix_t *orbitSizes = NULL;
static MtxFile_t *fileInp = NULL;
static MtxFile_t *fileOut = NULL;
static PTR rowInp = NULL;
static PTR rowOut = NULL;
static uint32_t degree;
static uint32_t nOrbits = 0;    // number of orbits (= input vector size)
static uint32_t nVectors = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////

static void readOrbits()
{
   MtxFile_t* orbitsFile = mfOpen(fileNameOrbits);
   orbits = imatRead(orbitsFile);
   orbitSizes = imatRead(orbitsFile);
   degree = orbits->noc;
   nOrbits = orbitSizes->noc;
   mfClose(orbitsFile);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void openFiles()
{
   // Vector input file.
   fileInp = mfOpen(fileNameInp);
   mfReadHeader(fileInp);
   if (mfObjectType(fileInp) != MTX_TYPE_MATRIX)
      mtxAbort(MTX_HERE,"%s: %s",fileNameInp,MTX_ERR_NOTMATRIX);
   if (fileInp->header[2] != nOrbits)
      mtxAbort(MTX_HERE,"%s and %s: %s",fileNameInp,fileNameOrbits,MTX_ERR_INCOMPAT);
   nVectors = fileInp->header[1];
   ffSetField(fileInp->header[0]);
   rowInp = ffAlloc(1, nOrbits);

    // Vector output file.
    fileOut = mfCreate(fileNameOut,ffOrder,nVectors,degree);
    rowOut = ffAlloc(1, degree);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void uncondense()
{
   for (uint32_t i = 0; i < nVectors; ++i) {
      mfReadRows(fileInp, rowInp, 1, nOrbits);

      ffMulRow(rowOut, FF_ZERO, degree);
      for (uint32_t k = 0; k < nOrbits; ++k) {
         int l;
         FEL f = ffExtract(rowInp, k);
         int count = orbitSizes->data[k];
         for (l = 0; count > 0 && l < degree; ++l) {
            if (orbits->data[l] == k) {
               ffInsert(rowOut, l, f);
               --count;
            }
         }
      }

      mfWriteRows(fileOut, rowOut, 1, degree);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
   app = appAlloc(&AppInfo,argc,argv);
   appGetArguments(app,3,3);
   fileNameInp = app->argV[0];
   fileNameOrbits = app->argV[1];
   fileNameOut = app->argV[2];
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void cleanup()
{
   mfClose(fileInp);
   mfClose(fileOut);
   imatFree(orbits);
   imatFree(orbitSizes);
   appFree(app);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   init(argc,argv);
   readOrbits();
   openFiles();
   uncondense();
   cleanup();
   return 0;
}




////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_zuk zuk - Uncondense Vectors

@section zuk_syntax Command Line
<pre>
zuk @em Options @em Vectors @em Orbits @em Result
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par @em Vectors
  Vectors to uncondense.
@par @em Orbits
  Orbit tables produced by @ref prog_zmo "zmo".
@par @em Result
  Uncondensed vectors.

@section zuk_inp Input Files
@par @em Vectors
  Vectors to uncondense.
@par @em Orbits
  Orbit tables produced by @ref prog_zmo "zmo".

@section zuk_out Output Files
@par @em Result
  Uncondensed vectors.

@section zuk_desc Description
This program reads a matrix which is assumed to be a condensed space of a permutation
representation whose orbits are in the file @em Orbits. The vectors in @em Vectors are
elongated so as to lie in the original permutation space and written out to 
the file @em Result.
@em Orbits must be an orbits file in the format defined by @ref prog_zmo "zmo".
Here is an example:
<pre>
         2 0 4
Space =  1 3 2
         2 0 2

Orbits = (1,2) (3,4,5,6) (7,8,9)

         2 2 0 0 0 0 4 4 4
Result = 1 1 3 3 3 3 2 2 2
         2 2 0 0 0 0 2 2 2
</pre>
*/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
