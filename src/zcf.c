////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Change field.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <stdlib.h>



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


static const char *iname, *oname;
static MtxFile_t* inputFile = NULL;
static MtxFile_t *outputFile = NULL;
static int inputFieldOrder;
static int outputFieldOrder;
//static int nor, noc;		    /* Parameters of input file */


static MtxApplicationInfo_t AppInfo = { 
"zcf", "Change Format", 
"SYNTAX\n"
"    zcf " MTX_COMMON_OPTIONS_SYNTAX " <Field> <Input> <Output>\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"\n"
"ARGUMENTS\n"
"    <Field> ................. Desired field order\n"
"    <Input> ................. Input file name\n"
"    <Output> ................ Output file name\n"
};

static MtxApplication_t *App = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////

static void checkFieldOrders()
{	
    if (outputFieldOrder < 2)
    	mtxAbort(MTX_HERE,"Invalid target field order %d",(int) outputFieldOrder);

    if (outputFieldOrder == inputFieldOrder) {
    	mtxAbort(MTX_HERE,"%s is already over GF(%d)",iname, (int) outputFieldOrder);
    }
       
    int *sfp = mtx_subfields;
    if (outputFieldOrder > inputFieldOrder) {
       ffSetField(outputFieldOrder);
       while (*sfp >= 2 && *sfp != inputFieldOrder) 
          ++sfp;
    } else {
       ffSetField(inputFieldOrder);
       while (*sfp >= 2 && *sfp != outputFieldOrder) 
          ++sfp;
    }
    if (*sfp < 2)
       mtxAbort(MTX_HERE,"Cannot change from GF(%d) to GF(%d)",inputFieldOrder,outputFieldOrder);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void convertPermutationToMatrix()
{
   ffSetField(outputFieldOrder);
   Perm_t *perm = permReadData(inputFile->file, inputFile->header);
   const uint32_t nor = perm->degree;
   PTR row = ffAlloc(1, nor);

   outputFile = mfCreate(oname, outputFieldOrder, nor, nor);
   const uint32_t *p = perm->data;
   for (uint32_t i = 0; i < nor; ++i)
   {	
      ffMulRow(row,FF_ZERO, nor);
      ffInsert(row,p[i],FF_ONE);
      mfWriteRows(outputFile,row,1, nor);
   }
   permFree(perm);
   sysFree(row);
   MESSAGE(0,("Converted to GF(%d)\n",outputFieldOrder));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void changeField()
{
   inputFieldOrder = inputFile->header[0];
   const uint32_t nor = inputFile->header[1];
   const uint32_t noc = inputFile->header[2];
   checkFieldOrders();

   // Allocate buffers
   static const size_t MAX_ROWS = 1000;       // chunk size
   FEL *bufTmp = NALLOC(FEL,MAX_ROWS * noc);  // MAX_ROWS * noc unpacked field elements
   ffSetField(inputFieldOrder);
   PTR bufInp = ffAlloc(MAX_ROWS, noc);       // MAX_ROWS rows
   ffSetField(outputFieldOrder);
   PTR bufOut = ffAlloc(MAX_ROWS, noc);       // MAX_ROWS rows

   outputFile = mfCreate(oname,outputFieldOrder,nor,noc);

   // Process data in chunks of MAX_ROWS rows.
   for (uint32_t rowsLeft = nor; rowsLeft > 0; )
   {
      // Read next chunk of rows into buffer.
      const uint32_t rowsRead = (rowsLeft <= MAX_ROWS) ? rowsLeft : MAX_ROWS;
      ffSetField(inputFieldOrder);
      mfReadRows(inputFile,bufInp,nor,noc);

      // Unpack input rows
      {
         PTR row = bufInp;
         FEL* bp = bufTmp;
         for (size_t r = rowsRead; r > 0; --r) {
            for (size_t c = 0; c < noc; ++c) {
               *bp++ = ffExtract(row, c);
            }
            ffStepPtr(&row, noc);
         }
      }

      // Convert to target field.
      FEL* const bufEnd = bufTmp + (size_t) rowsRead * noc;
      if (inputFieldOrder < outputFieldOrder) {
         ffSetField(outputFieldOrder);
         for (FEL * bp = bufTmp; bp < bufEnd; ++bp)
            *bp = ffEmbed(*bp, inputFieldOrder);
      } else {
         ffSetField(inputFieldOrder);
         for (FEL * bp = bufTmp; bp < bufEnd; ++bp)
            *bp = ffRestrict(*bp, outputFieldOrder);
      }

      // Pack into output buffer.
      ffSetField(outputFieldOrder);
      {
         PTR row = bufOut;
         FEL* bp = bufTmp;
         for (size_t r = rowsRead; r > 0; --r) {
            for (size_t c = 0; c < noc; ++c) {
               ffInsert(row, c, *bp++);
            }
            ffStepPtr(&row, noc);
         }
      }

      // Write output rows
      mfWriteRows(outputFile, bufOut, rowsRead, noc);

      rowsLeft -= rowsRead;
   }

   sysFree(bufOut);
   sysFree(bufTmp);
   sysFree(bufInp);

   if (inputFieldOrder < outputFieldOrder)
      MESSAGE(0,("Embedded into GF(%d)\n",outputFieldOrder));
   else
      MESSAGE(0,("Restricted to GF(%d)\n",outputFieldOrder));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    appGetArguments(App,3,3);
    outputFieldOrder = atol(App->argV[0]);
    iname = App->argV[1];
    oname = App->argV[2];
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void cleanUp()
{
    if (App != NULL)
	appFree(App);
    if (inputFile != NULL)
	mfClose(inputFile);
    if (outputFile != NULL)
	mfClose(outputFile);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)

{
    init(argc,argv);

    inputFile = mfOpen(iname);
    mfReadHeader(inputFile);
    uint32_t objType = mfObjectType(inputFile);
    if (objType == MTX_TYPE_MATRIX)
       changeField();
    else if (objType == MTX_TYPE_PERMUTATION)
       convertPermutationToMatrix();
    else
       mtxAbort(MTX_HERE, "%s: unsupported object type", iname);

    cleanUp();
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_zcf zcf - Change Field

@section zcf_syntax Command Line
<pre>
zcf @em Options @em q @em Input @em Output
</pre>

@par @em Options
Standard options, see @ref prog_stdopts

@par @em q
New field order.

@par @em Input
Input file.

@par @em Output
Output file.

@section zcf_inp Input Files

@par @em Input
Input file.

@section zcf_out Output Files

@par @em Output
Output file.

@section zcf_desc Description

This program converts between various data types. Currently
there are two kinds of conversions available:
- If @em Input is a matrix, the field is changed to GF(@em q).
  The current field of the input file must be a subfield
  or a superfield of GF(@em q). In the latter case, all matrix
  entries must be in GF(@em q).
- If @em Input is a permutation of degree n, it is converted
  into the corresponding n times n permutation matrix over
  GF(@em q).

@section zcf_impl Implementation Details
For matrices, the conversion is done in two steps. First,
all entries of the matrix are converted to integers.
Then, they are mapped to the new field and reassembled into rows.
The result is written outputFile row by row.

In case of permutations the output matrix is generated row by
row by inserting ones at the positions specified by the
permutation.

If the input is a matrix, the whole matrix must fit into memory.
Additionally, the program needs n⋅m⋅s bytes of memory, where m and n
are the dimensions of the input matrix and s=1 for the small
arithmetic version and s=2 for the big version.
In case of permutations, the input permutation and one row
of the output file must fit into memory.
*/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
