////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - print a matrix or permutaion.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

static FILE* textFile = NULL;
static MtxFile_t* binaryFile = NULL;
static int Gap = 0;
static int Summary = 0;

static MtxApplicationInfo_t AppInfo = {
   "zpr", "print Permutations Or Matrices",
   "SYNTAX\n"
   "    zpr [-G] [-s] <Binfile> [<Textfile>]\n"
   "\n"
   "OPTIONS\n"
   "    -G   GAP output\n"
   "    -s   print summary only\n"
   "\n"
   "FILES\n"
   "    <Binfile>   i  A matrix or permutation in binary format\n"
   "    <Textfile>  i  The output in text format (default: stdout)\n"
};

static MtxApplication_t* App = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////

static unsigned width = 0;
static unsigned maxWidth = 80;

static void printNewLine()
{
   if (width > 0) {
      fputc('\n', textFile);
      width = 0;
   }
}

MTX_PRINTF_ATTRIBUTE(1,2)
static void print(const char *fmt, ...)
{
   char tmp[200];
   va_list args;
   va_start(args, fmt);
   const int n = vsnprintf(tmp, sizeof(tmp), fmt, args);
   va_end(args);
   MTX_ASSERT(n >= 0);
   if (width + n > maxWidth) {
      printNewLine();
   }
   fwrite(tmp, 1, n, textFile);
   width += n;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void printString(char* c)
{
   static int pos = 0;
   const int l = strlen(c);

   if (l + pos >= 78) {
      fprintf(textFile, "\n");
      pos = 0;
   }
   fprintf(textFile, "%s", c);
   for (; *c != 0; ++c) {
      if (*c == '\n') {
         pos = 0;
      }
      else {
         ++pos;
      }
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void print32u(uint32_t value)
{
   char tmp[50];
   snprintf(tmp, sizeof(tmp), "%lu", (unsigned long) value);
   printString(tmp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void print32s(int32_t value)
{
   char tmp[50];
   snprintf(tmp, sizeof(tmp), "%ld", (unsigned long) value);
   printString(tmp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void printGapMatrix()
{
   const uint32_t field = binaryFile->header[0];
   const uint32_t nor = binaryFile->header[1];
   const uint32_t noc = binaryFile->header[2];

   ffSetField(field);
   const int isprimefield = (ffChar == ffOrder);
   PTR m1 = ffAlloc(1, noc);
   printString("MeatAxe.Matrix := [\n");
   for (uint32_t row = 1; row <= nor; ++row) {
      mfReadRows(binaryFile, m1, 1, noc);
      int cnt = 0;
      fprintf(textFile, "[");
      for (uint32_t col = 0; col < noc; ++col) {
         if (cnt > 75) {
            fprintf(textFile, "\n ");
            cnt = 0;
         }
         const FEL f1 = ffExtract(m1, col);
         if (isprimefield) {
            FEL f2 = FF_ZERO;
            unsigned long k = 0;
            while (f2 != f1) {
               f2 = ffAdd(f2, ffGen);
               ++k;
            }
            cnt += fprintf(textFile, "%lu", k);
         } else {
            if (f1 == FF_ZERO) {
               cnt += fprintf(textFile, "0*Z(%lu)", (unsigned long)field);
            } else {
               FEL f2 = ffGen;
               unsigned long k = 1;
               while (f2 != f1) {
                  f2 = ffMul(f2, ffGen);
                  ++k;
               }
               cnt += fprintf(textFile, "Z(%lu)^%lu", (unsigned long)field, k);
            }
         }
         if (col < noc - 1) {
            fprintf(textFile, ",");
            ++cnt;
         }
      }
      fprintf(textFile, "]");
      if (row < nor) {
         fprintf(textFile, ",");
      }
      fprintf(textFile, "\n");
   }
   fprintf(textFile, "]");
   if (isprimefield) {
      fprintf(textFile, "*Z(%lu)", (unsigned long)field);
   }
   fprintf(textFile, ";\n");

   sysFree(m1);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void printGapIntegerMatrix()
{
   int32_t* row;
   const uint32_t nor = binaryFile->header[1];
   const uint32_t noc = binaryFile->header[2];

   row = NALLOC(int32_t, noc);
   printString("MeatAxe.Matrix := [\n");
   for (uint32_t r = 0; r < nor; ++r) {
      mfRead32(binaryFile, row, noc);
      fprintf(textFile, "[");
      for (uint32_t c = 0; c < noc; ++c) {
         int32_t k = row[c];
         print32s(k);
         if (c < noc - 1) {
            fprintf(textFile, ",");
         }
      }
      fprintf(textFile, "]");
      if (r + 1 < nor) {
         fprintf(textFile, ",");
      }
      fprintf(textFile, "\n");
   }
   fprintf(textFile, "];\n");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void printGapPermutation()
{
   const uint32_t numberOfPermutations = binaryFile->header[2];
   const uint32_t degree = binaryFile->header[1];
   uint32_t* perm = NALLOC(uint32_t, degree);

   printString("MeatAxe.Perms := [\n");
   for (uint32_t pos = 1; pos <= numberOfPermutations; ++pos) {
      mfRead32(binaryFile, perm, degree);

      printString("    PermList([");
      for (uint32_t i = 0; i < degree; ++i) {
         if (i > 0) {
            printString(",");
         }
         print32u(perm[i] + 1);
      }
      printString("])");
      if (pos < numberOfPermutations) { printString(",");}
      printString("\n");
   }
   printString("];\n");
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void printGapFormat()
{
   const uint32_t objectType = mfObjectType(binaryFile);
   if (objectType == MTX_TYPE_PERMUTATION) {
      printGapPermutation();
   } else if (objectType == MTX_TYPE_INTMATRIX)
   {
      printGapIntegerMatrix();
   } else if (objectType == MTX_TYPE_MATRIX)
   {
      printGapMatrix();
   } else {
      mtxAbort(MTX_HERE, "Cannot print type 0x%lu in GAP format", (unsigned long)objectType);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void printMatrix()
{
   int width, marksPerLine;
   const uint32_t q = binaryFile->header[0];
   const uint32_t nor = binaryFile->header[1];
   const uint32_t noc = binaryFile->header[2];
   ffSetField(q);
   PTR m1 = ffAlloc(1, noc);

   if (q < 10) {
      width = 1;
      marksPerLine = 80;
   }
   else if (q < 100) {
      width = 3;
      marksPerLine = 25;
   }
   else if (q < 1000) {
      width = 4;
      marksPerLine = 20;
   }
   else if (q < 10000) {
      width = 5;
      marksPerLine = 15;
   }
   else {
      width = 6;
      marksPerLine = 12;
   }

   fprintf(textFile, "matrix field=%lu rows=%lu cols=%lu\n",
           (unsigned long)q, (unsigned long)nor, (unsigned long)noc);
   for (uint32_t row = 0; row < nor; ++row) {
      ffReadRows(binaryFile->file, m1, 1, noc);
      for (int c = 0; c < noc; ++c) {
         const int value = ffToInt(ffExtract(m1, c));
         fprintf(textFile, "%*d", width, value);
         if ((c + 1) % marksPerLine == 0 || c + 1 == noc) {
            fprintf(textFile, "\n");
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void printPolynomial()
{
   Poly_t* p = polReadData(binaryFile->file, binaryFile->header);
   print("polynomial field=%lu degree=%ld", (unsigned long)p->field, (long)p->degree);
   printNewLine();
   for (int32_t i = 0; i <= p->degree; ++i) {
      print("%s%lu", i > 0 ? " " : "", (unsigned long) ffToInt(p->data[i]));
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void printPermutation()
{
   Perm_t* perm = permReadData(binaryFile->file, binaryFile->header);
   print("permutation degree=%lu", (unsigned long) perm->degree);
   printNewLine();
   for (uint32_t i = 0; i < perm->degree; ++i) {
      print("%s%lu", i > 0 ? " " : "", (unsigned long)perm->data[i] + 1);
   }
   permFree(perm);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void printIntegerMatrix()
{
   const uint32_t nor = binaryFile->header[1];
   const uint32_t noc = binaryFile->header[2];
   int32_t* row = NALLOC(int32_t, noc);
   print("integer-matrix rows=%lu cols=%lu", (unsigned long)nor, (unsigned long)noc);
   printNewLine();
   for (uint32_t i = 0; i < nor; ++i) {
      mfRead32(binaryFile, row, noc);
      for (uint32_t k = 0; k < noc; ++k) {
         print("%s%ld", k > 0 ? " " : "", (long) row[k]);
      }
      printNewLine();
   }
   sysFree(row);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void printStandardFormat()
{
   const uint32_t objectType = mfObjectType(binaryFile);
   if (objectType == MTX_TYPE_MATRIX) {
      printMatrix();
   }
   else if (objectType == MTX_TYPE_PERMUTATION) {
      printPermutation();
   }
   else if (objectType == MTX_TYPE_POLYNOMIAL) {
      printPolynomial();
   }
   else if (objectType == MTX_TYPE_INTMATRIX) {
      printIntegerMatrix();
   }
   else {
      mtxAbort(MTX_HERE, "Cannot print type 0x%lx in Mtx format", (unsigned long)objectType);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void printPermutationSummary()
{
   const uint32_t degree = binaryFile->header[1];
   const uint32_t nPerms = binaryFile->header[2];
   if (Gap) {
      printf("MeatAxe.PermutationCount:=%ld;\n", (long)nPerms);
      printf("MeatAxe.PermutationDegree:=%ld;\n", (long)degree);
   } else {
      printf("%ld Permutation%s of degree %ld\n",
             (long)nPerms, nPerms == 1 ? "" : "s", (long)degree);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void printMatrixSummary()
{
   const uint32_t field = binaryFile->header[0];
   const uint32_t nor = binaryFile->header[1];
   const uint32_t noc = binaryFile->header[2];
   if (Gap) {
      printf("MeatAxe.MatrixRows:=%ld;\n", (long)nor);
      printf("MeatAxe.MatrixCols:=%ld;\n", (long)noc);
      printf("MeatAxe.MatrixField:=%ld;\n", (long)field);
   } else {
      printf("%ld x %ld matrix over GF(%ld)\n", (long)nor, (long)noc, (long)field);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void printPolySummary()
{
   const uint32_t field = binaryFile->header[1];
   const int32_t degree = (int32_t) binaryFile->header[2];
   if (Gap) {
      printf("MeatAxe.PolynomialField=%lu\n", (unsigned long)field);
      printf("MeatAxe.PolynomialDegree=%ld\n", (long)degree);
   } else {
      printf("Polynomial of degree %ld over GF(%lu)\n", (long)degree, (unsigned long)field);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void printImatSummary()
{
   const uint32_t nor = binaryFile->header[1];
   const uint32_t noc = binaryFile->header[2];
   if (Gap) {
      printf("MeatAxe.IntegerMatrixRows:=%lu;\n", (unsigned long)nor);
      printf("MeatAxe.IntegerMatrixCols:=%lu;\n", (unsigned long)noc);
   }
   else {
      printf("%lu x %lu integer matrix\n", (unsigned long)nor, (unsigned long)noc);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void printSummary()
{
   const uint32_t objectType = mfObjectType(binaryFile);
   size_t objectSize = 0;
   if (objectType == MTX_TYPE_PERMUTATION) {
      printPermutationSummary();
      objectSize = sizeof(uint32_t) * binaryFile->header[1] * binaryFile->header[2];
   }
   else if (objectType == MTX_TYPE_MATRIX) {
      printMatrixSummary();
      ffSetField(binaryFile->header[0]);
      objectSize = ffRowSizeUsed(binaryFile->header[2]) * binaryFile->header[1];
   }
   else if (objectType == MTX_TYPE_POLYNOMIAL) {
      printPolySummary();
      ffSetField(binaryFile->header[1]);
      const int32_t degree = (int32_t)binaryFile->header[2];
      MTX_ASSERT(degree >= -1);
      objectSize = ffRowSizeUsed(degree + 1);
   }
   else if (objectType == MTX_TYPE_INTMATRIX) {
      printImatSummary();
      objectSize = sizeof(int32_t) * binaryFile->header[1] * binaryFile->header[2];
   }
   else {
      mtxAbort(MTX_HERE, "Unsupported/invalid file header (0x%lx,0x%lx,0x%lx)",
               (unsigned long)binaryFile->header[0],
               (unsigned long)binaryFile->header[1],
               (unsigned long)binaryFile->header[2]);
   }

   mfSkip(binaryFile, objectSize);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char** argv)
{
   App = appAlloc(&AppInfo, argc, argv);
   Gap = appGetOption(App, "-G --gap");
   Summary = appGetOption(App, "-s --summary");
   if (Gap) {
      MtxMessageLevel = -100;   /* Suppress messages in GAP mode */
   }
   appGetArguments(App, 1, 2);
   binaryFile = mfOpen(App->argV[0]);
   textFile = (App->argC >= 2) ? sysFopen(App->argV[1], "w") : stdout;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
   init(argc, argv);
   while (mfTryReadHeader(binaryFile)) {
      if (Summary) {
         printSummary();
      }
      else if (Gap) {
         printGapFormat();
      }
      else {
         printStandardFormat();
      }
      printNewLine();
   }
   mfClose(binaryFile);
   fclose(textFile);
   return EXIT_OK;
}

// @page prog_zpr zpr - Print Matrices and Permutations
// @see @ref prog_zcv
//
// @section zpr_syntax Command Line
// <pre>
// zpr [@em Options] [-Gs] @em DataFile [@em TextFile]
// </pre>
//
// @par @em Options
// Standard options, see @ref prog_stdopts
// @par -G, --gap
// Output in GAP format.
// @par -s, --summary
// Show headers only.
// @par @em DataFile
// Input file (binary)
// @par @em TextFile
// Output file (text)
//
// @section zpr_inp Input Files
// @par @em DataFile
// Input file (binary)
//
// @section zpr_out Output Files
// @par @em TextFile
// Output file (text)
//
// @section zpr_desc Description
// This program prints the contents of a MeatAxe data file in readable
// format. The text produced by @b zpr can be converted into binary format by
// the @ref prog_zcv "zcv" program.
//
// If there is only one argument on the command line, @b zpr writes to stdout.
// A second argument, if present, is taken as the output file name.
//
// To find out the contents of a MeatAxe file, use the -s option. To generate
// output readble by GAP, use -G. Both options can be combined.

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
