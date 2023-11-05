////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Sum of matrices.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

/* ------------------------------------------------------------------
   Variables
   ------------------------------------------------------------------ */

static MtxApplicationInfo_t AppInfo = {
   "zad", "Add or Subtract Matrices",
   "SYNTAX\n"
   "    zad " MTX_COMMON_OPTIONS_SYNTAX " [-]<Mat> [-]<Mat> ... <Result>\n"
   "\n"
   "ARGUMENTS\n"
   "    <Mat> ................... Input file: Matrix to add (-<Mat> subtracts)\n"
   "    <Result> ................ Output file: Sum\n"
   "\n"
   "OPTIONS\n"
   MTX_COMMON_OPTIONS_DESCRIPTION
   "\n"
   "FILES\n"
   "    <Mat> ................... I Input matrix\n"
   "    <Result> ................ O Sum of the input matrices\n"
};

static MtxApplication_t* App = NULL;

#define MTX_MAX_INPUT 20
int NInput;                         /* Number of input matrices */
MtxFile_t* Input[MTX_MAX_INPUT] = { 0 };
int Subtract[MTX_MAX_INPUT];
uint32_t Field = 0;
uint32_t Nor = 0;
uint32_t Noc = 0;
PTR Buf1, Buf2;                     /* Working buffer */
MtxFile_t* Output = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char** argv)
{
   App = appAlloc(&AppInfo, argc, argv);
   appGetArguments(App, 3, MTX_MAX_INPUT + 1);
   NInput = App->argC - 1;

   // Open the input files
   for (int i = 0; i < NInput; ++i) {
      const char* file_name = App->argV[i];
      Subtract[i] = *file_name == '-';
      if (*file_name == '-' || *file_name == '+') {
         ++file_name;
      }
      Input[i] = mfOpen(file_name);
      mfReadHeader(Input[i]);

      if (mfObjectType(Input[i]) != MTX_TYPE_MATRIX) {
         mtxAbort(MTX_HERE, "%s: %s (type=0x%lx)",
                  Input[i]->name,
                  MTX_ERR_NOTMATRIX,
                  (long)Input[i]->header[0]);
      }
      if (i == 0) {
         Field = Input[i]->header[0];
         Nor = Input[i]->header[1];
         Noc = Input[i]->header[2];
      }
      else if (Input[i]->header[0] != Field
               || Input[i]->header[1] != Nor || Input[i]->header[2] != Noc) {
         mtxAbort(MTX_HERE, "%s and %s: %s", Input[0]->name, Input[i]->name, MTX_ERR_INCOMPAT);
      }
   }

   // Open the output file.
   Output = mfCreate(App->argV[App->argC - 1], Field, Nor, Noc);

   // Allocate workspace
   ffSetField(Field);
   Buf1 = ffAlloc(1, Noc);
   Buf2 = ffAlloc(1, Noc);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void cleanUp()
{
   // Close all files.
   for (int i = 0; i < NInput; ++i) {
      mfClose(Input[i]);
   }
   mfClose(Output);

   // Free workspace.
   sysFree(Buf1);
   sysFree(Buf2);
   appFree(App);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void addMatrices()
{
   const FEL MinusOne = ffNeg(FF_ONE);

   for (uint32_t i = Nor; i > 0; --i) {
      mfReadRows(Input[0], Buf1, 1, Noc);
      if (Subtract[0]) {
         ffMulRow(Buf1, MinusOne, Noc);
      }
      for (int k = 1; k < NInput; ++k) {
         mfReadRows(Input[k], Buf2, 1, Noc);
         if (Subtract[k]) {
            ffAddMulRow(Buf1, Buf2, MinusOne, Noc);
         }
         else {
            ffAddRow(Buf1, Buf2, Noc);
         }
      }

      mfWriteRows(Output, Buf1, 1, Noc);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
   init(argc, argv);
   addMatrices();
   cleanUp();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_zad zad - Add Matrices

@section zad_syntax Command Line
<pre>
zad [@em Options] [-]@em Mat [-]@em Mat ... @em Result
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par @em Mat
  Input matrix
@par @em Result
  Result matrix

@section zad_inp Input Files
@par @em Mat
  Input matrix

@section zad_out Output Files
@par @em Result
  Result matrix

@section zad_desc Description

This program reads two or more input matrices, calculates their sum or difference and
writes the result to a file.  The input matrices must be compatible, i.e., they must be
over the same field and have the same dimensions. @b zad is designed to work with very
large matrices without running out of memory. Only two rows are allocated as working 
memory.

By default, all input matrices are added. For example,
<pre>
zad A B C
</pre>
calculates the sum of A and B, and writes the result to C.
If a file name is preceeded by a minus sign, this matrix is subtracted.
For example,
<pre>
zad A -B -C D
</pre>
calculates $D=A-B-C$. To subtract the first matrix, you must insert an
extra "--" before the file names. Otherwise the first argument
would be interpreted as a program option. For example, to calculate $C=-A-B$, use
<pre>
zad -- -A -B C
</pre>
If a file name starts with "-", preceed the file name by "+" to add, or
by "-" to subtract the matrix. If, for example, the second input file
is "-B", use the following syntax:
<pre>
zad A +-B C
zad A --B C
</pre>

*/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
