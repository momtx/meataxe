////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Print a matrix or permutaion.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static FILE* textFile = NULL;   /* Output file */
static const char* binaryFileName = NULL;
static MtxFile_t* binaryFile = NULL;
static int Gap = 0;             /* -g (GAP mode) */
static int Summary = 0;         /* -s (summary) */

static MtxApplicationInfo_t AppInfo = {
   "zpr", "Print Permutations Or Matrices",
   "SYNTAX\n"
   "    zpr [-G] [-s] <Binfile> [<Textfile>]\n"
   "\n"
   "OPTIONS\n"
   "    -G   GAP output\n"
   "    -s   Print summary only\n"
   "\n"
   "FILES\n"
   "    <Binfile>   i  A matrix or permutation in binary format\n"
   "    <Textfile>  i  The output in text format (default: stdout)\n"
};

static MtxApplication_t* App = NULL;

/* ------------------------------------------------------------------
   PrString(), PrLong() - Prettty printer
   ------------------------------------------------------------------ */

static void PrString(char* c)
{
   static int pos = 0;
   int l = strlen(c);

   if (l + pos >= 78) {
      fprintf(textFile,"\n");
      pos = 0;
   }
   fprintf(textFile,"%s",c);
   for (; *c != 0; ++c) {
      if (*c == '\n') { pos = 0;} else { ++pos;}}
}


static void PrLong(long l)
{
   char tmp[50];
   snprintf(tmp, sizeof(tmp), "%ld",l);
   PrString(tmp);
}


/* ------------------------------------------------------------------
   prmatrix() - Print a matrix in standard format.
   ------------------------------------------------------------------ */

static void prmatrix()
{
   int md, mx;
   const uint32_t q = binaryFile->header[0];
   const uint32_t nor = binaryFile->header[1];
   const uint32_t noc = binaryFile->header[2];
   ffSetField(q);
   PTR m1 = ffAlloc(1, noc);

   if (q < 10) {md = 1; mx = 80;} else if (q < 100) {md = 3; mx = 25;} else if (q < 1000) {
      md = 4;
      mx = 20;
   } else if (q < 10000) {md = 5; mx = 15;} else {md = 6; mx = 12;}

   fprintf(textFile,"matrix field=%ld rows=%ld cols=%ld\n", (long)q,(long)nor,(long)noc);
   for (int row = 1; row <= nor; ++row) {
      ffReadRows(binaryFile->File, m1, 1, noc);
      int iv = 1;
      for (int j1 = 0; j1 < noc; ++j1) {
         FEL f1 = ffExtract(m1,j1);
         switch (md) {
            case 1: fprintf(textFile,"%1d",ffToInt(f1));
               break;
            case 2: fprintf(textFile,"%2d",ffToInt(f1));
               break;
            case 3: fprintf(textFile,"%3d",ffToInt(f1));
               break;
            case 4: fprintf(textFile,"%4d",ffToInt(f1));
               break;
            case 5: fprintf(textFile,"%5d",ffToInt(f1));
               break;
            case 6: fprintf(textFile,"%6d",ffToInt(f1));
               break;
         }
         if (iv++ >= mx) {
            fprintf(textFile,"\n");
            iv = 1;
         }
      }
      if (iv > 1) {
         fprintf(textFile,"\n");
      }
   }
}


/* ------------------------------------------------------------------
   prgapmat() - Print a matrix in GAP format.
   ------------------------------------------------------------------ */

static int nDigits(int k)
{
   if (k >= 100000) return 6;
   if (k >= 10000) return 5;
   if (k >= 1000) return 4;
   if (k >= 100) return 3;
   if (k >= 10) return 2;
   return 1;
}

static void prgapmat()
{
   const uint32_t field = binaryFile->header[0];
   const uint32_t nor = binaryFile->header[1];
   const uint32_t noc = binaryFile->header[2];

   ffSetField(field);
   const int isprimefield = (ffChar == ffOrder);
   PTR m1 = ffAlloc(1, noc);
   PrString("MeatAxe.Matrix := [\n");
   for (uint32_t loop1 = 1; loop1 <= nor; ++loop1) {
      mfReadRows(binaryFile, m1, 1, noc);
      int cnt = 0;
      fprintf(textFile,"[");
      for (uint32_t col = 0; col < noc; ++col) {
         if (cnt > 75) {
            fprintf(textFile,"\n ");
            cnt = 0;
         }
         const FEL f1 = ffExtract(m1,col);
         if (isprimefield) {
            FEL f2 = FF_ZERO;
            long k = 0;
            while (f2 != f1) {
               f2 = ffAdd(f2,ffGen);
               ++k;
            }
            fprintf(textFile,"%ld",k);
            cnt += nDigits(k);
         } else {
            if (f1 == FF_ZERO) {
               fprintf(textFile,"0*Z(%ld)",(long)field);
               cnt += 5;
               cnt += nDigits(field);
            } else {
               FEL f2 = ffGen;
               long k = 1;
               while (f2 != f1) {
                  f2 = ffMul(f2,ffGen);
                  ++k;
               }
               fprintf(textFile,"Z(%ld)^%ld",(long)field,k);
               cnt += 4 + nDigits(field) + nDigits(k);
            }
         }
         if (col < noc - 1) {
            fprintf(textFile,",");
            ++cnt;
         }
      }
      fprintf(textFile,"]");
      if (loop1 < nor) {
         fprintf(textFile,",");
      }
      fprintf(textFile,"\n");
   }
   fprintf(textFile,"]");
   if (isprimefield) {
      fprintf(textFile,"*Z(%ld)",(long)field);
   }
   fprintf(textFile,";\n");

   sysFree(m1);
}


/* ------------------------------------------------------------------
   prgapimat() - Print an integer matrix in GAP format.
   ------------------------------------------------------------------ */

static void prgapimat()
{
   uint32_t* row;
   long loop1, j1;
   int cnt;
   const uint32_t nor = binaryFile->header[1];
   const uint32_t noc = binaryFile->header[2];

   row = NALLOC(uint32_t,noc);
   PrString("MeatAxe.Matrix := [\n");
   for (loop1 = 1; loop1 <= nor; ++loop1) {
      mfRead32(binaryFile,row,noc);
      cnt = 0;
      fprintf(textFile,"[");
      for (j1 = 0; j1 < noc; ++j1) {
         long k = row[j1];
         if (cnt > 75) {
            fprintf(textFile,"\n ");
            cnt = 0;
         }
         fprintf(textFile,"%ld",k);
         cnt += k > 9999?5:k > 999?4:k > 99?3:k > 9?2:1;
         if (j1 < noc - 1) {
            fprintf(textFile,",");
            ++cnt;
         }
      }
      fprintf(textFile,"]");
      if (loop1 < nor) {
         fprintf(textFile,",");
      }
      fprintf(textFile,"\n");
   }
   fprintf(textFile,"];\n");
}


/* ------------------------------------------------------------------
   prgapperm() - Print a permutation in GAP format.
   ------------------------------------------------------------------ */

static void prgapperm()
{
   long i, pos;
   uint32_t* perm;

   const size_t numberOfPermutations = binaryFile->header[2];
   const size_t degree = binaryFile->header[1];

   perm = NALLOC(uint32_t, degree);
   if (perm == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate work space");
   }

   PrString("MeatAxe.Perms := [\n");
   for (pos = 1; pos <= numberOfPermutations; ++pos) {
      // Read the next permutation
      mfRead32(binaryFile,perm,degree);

      PrString("    PermList([");
      for (i = 0; i < degree; ++i) {
         if (i > 0) {
            PrString(",");
         }
         PrLong(perm[i] + 1);
      }
      PrString("])");
      if (pos < numberOfPermutations) { PrString(",");}
      PrString("\n");
   }
   PrString("];\n");
}


/* ------------------------------------------------------------------
   prgap() - Print a matrix or permutation in GAP format.
   ------------------------------------------------------------------ */

static void prgap()
{
   if (binaryFile->header[0] == MTX_TYPE_PERMUTATION) {
      prgapperm();
   } else if (binaryFile->header[0] == MTX_TYPE_INTMATRIX) {
      prgapimat();
   } else if (binaryFile->header[0] >= 2) {
      prgapmat();
   } else {
      mtxAbort(MTX_HERE,"Cannot print type 0x%lu in GAP format",
               (unsigned long)binaryFile->header[0]);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void printPolynomial()
{
   Poly_t* p = polReadData(binaryFile->File, binaryFile->header);
   fprintf(textFile,"polynomial field=%ld degree=%ld\n",
           (unsigned long)p->Field, (unsigned long)p->Degree);
   for (uint32_t i = 0; i <= p->Degree; ++i) {
      PrLong(ffToInt(p->Data[i]));
      PrString(" ");
   }
   PrString("\n");
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void printPermutation()
{
   Perm_t* perm = permReadData(binaryFile->File, binaryFile->header);
   fprintf(textFile,"permutation degree=%d\n",perm->Degree);
   for (uint32_t i = 0; i < perm->Degree; ++i) {
      PrLong(perm->Data[i] + 1);
      PrString(" ");
   }
   PrString("\n");
   permFree(perm);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void printIntegerMatrix()
{
   const uint32_t nor = binaryFile->header[1];
   const uint32_t noc = binaryFile->header[2];
   int32_t* row = NALLOC(int32_t, noc);

   fprintf(textFile,"integer-matrix rows=%lu cols=%lu\n", (unsigned long)nor, (unsigned long)noc);
   for (uint32_t i = 0; i < nor; ++i) {
      mfRead32(binaryFile,row,noc);
      for (uint32_t k = 0; k < noc; ++k) {
         if (k > 0) PrString(" ");
         PrLong(row[k]);
      }
      PrString("\n");
   }
   sysFree(row);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void printStandardFormat()
{
   const uint32_t objectType = mfObjectType(binaryFile);
   if (objectType == MTX_TYPE_MATRIX) {
      prmatrix();
   } else if (objectType == MTX_TYPE_PERMUTATION) {
      printPermutation();
   } else if (objectType == MTX_TYPE_POLYNOMIAL) {
      printPolynomial();
   } else if (objectType == MTX_TYPE_INTMATRIX) {
      printIntegerMatrix();
   } else {
      mtxAbort(MTX_HERE,"Cannot print type 0x%lx in Mtx format", (unsigned long)objectType);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void PrintPermutationSummary()
{
   const uint32_t degree = binaryFile->header[1];
   const uint32_t nPerms = binaryFile->header[2];
   if (Gap) {
      printf("MeatAxe.PermutationCount:=%ld;\n",(long)nPerms);
      printf("MeatAxe.PermutationDegree:=%ld;\n",(long)degree);
   } else {
      printf("%ld Permutation%s of degree %ld\n",
            (long)nPerms, nPerms == 1 ? "" : "s", (long)degree);
   }
}


static void PrintMatrixSummary()
{
   const uint32_t field = binaryFile->header[0];
   const uint32_t nor = binaryFile->header[1];
   const uint32_t noc = binaryFile->header[2];
   if (Gap) {
      printf("MeatAxe.MatrixRows:=%ld;\n",(long)nor);
      printf("MeatAxe.MatrixCols:=%ld;\n",(long)noc);
      printf("MeatAxe.MatrixField:=%ld;\n",(long)field);
   } else {
      printf("%ld x %ld matrix over GF(%ld)\n", (long)nor, (long)noc, (long)field);
   }
}


static void PrintPolySummary()
{
   const uint32_t field = binaryFile->header[1];
   const uint32_t degree = binaryFile->header[2];
   if (Gap) {
      printf("MeatAxe.PolynomialField=%ld\n",(long)field);
      printf("MeatAxe.PolynomialDegree=%ld\n",(long)degree);
   } else {
      printf("Polynomial of degree %ld over GF(%ld)\n",(long)degree,(long)field);
   }
}


static void PrintImatSummary()
{
   const uint32_t nor = binaryFile->header[1];
   const uint32_t noc = binaryFile->header[2];
   if (Gap) {
      printf("MeatAxe.IntegerMatrixRows:=%ld;\n",(long)nor);
      printf("MeatAxe.IntegerMatrixCols:=%ld;\n",(long)noc);
   } else {
      printf("%ld x %ld integer matrix\n",(long)nor,(long)noc);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void PrintSummary()
{
   const uint32_t objectType = mfObjectType(binaryFile);
   size_t nSkip = 0;
   if (objectType == MTX_TYPE_PERMUTATION) {
      PrintPermutationSummary();
      nSkip = sizeof(uint32_t) * binaryFile->header[1] * binaryFile->header[2];
   } else if (objectType == MTX_TYPE_MATRIX) {
      PrintMatrixSummary();
      ffSetField(binaryFile->header[0]);
      nSkip = ffRowSizeUsed(binaryFile->header[2]) * binaryFile->header[1];
   } else if (objectType == MTX_TYPE_POLYNOMIAL) {
      PrintPolySummary();
      ffSetField(binaryFile->header[1]);
      nSkip = ffRowSizeUsed((size_t)binaryFile->header[2] + 1);
   } else if (objectType == MTX_TYPE_INTMATRIX) {
      PrintImatSummary();
      nSkip = sizeof(int32_t) * binaryFile->header[1] * binaryFile->header[2];
   } else {
      mtxAbort(MTX_HERE, "Unsupported/invalid file header (0x%lx,0x%lx,0x%lx)",
               (unsigned long)binaryFile->header[0],
               (unsigned long)binaryFile->header[1],
               (unsigned long)binaryFile->header[2]);
   }

   mfSkip(binaryFile, nSkip);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char** argv)
{
   // Process command line
   App = appAlloc(&AppInfo,argc,argv);
   Gap = appGetOption(App,"-G --gap");
   Summary = appGetOption(App,"-s --summary");
   if (Gap) {
      MtxMessageLevel = -100;   /* Suppress messages in GAP mode */
   }
   appGetArguments(App,1,2);

   binaryFileName = App->ArgV[0];
   binaryFile = mfOpen(binaryFileName);
   if (App->ArgC >= 2) {
      textFile = sysFopen(App->ArgV[1],"w");
   } else {
      textFile = stdout;
   }
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char** argv)
{
   init(argc,argv);
   while (mfTryReadHeader(binaryFile)) {
      if (Summary) {
         PrintSummary();
      } else if (Gap) {
         prgap();
      } else {
         printStandardFormat();
      }
   }
   mfClose(binaryFile);
   fclose(textFile);
   return EXIT_OK;
}


/**
   @page prog_zpr zpr - Print Matrices and Permutations
   @see @ref prog_zcv

   @section zpr_syntax Command Line
   <pre>
   zpr [@em Options] [-Gs] @em DataFile [@em TextFile]
   </pre>

   @par @em Options
   Standard options, see @ref prog_stdopts
   @par -G, --gap
   Output in GAP format.
   @par -s, --summary
   Show headers only.
   @par @em DataFile
   Input file (binary)
   @par @em TextFile
   Output file (text)

   @section zpr_inp Input Files
   @par @em DataFile
   Input file (binary)

   @section zpr_out Output Files
   @par @em TextFile
   Output file (text)

   @section zpr_desc Description
   This program prints the contents of a MeatAxe data file in readable
   format. The text produced by @b zpr can be converted into binary format by
   the @ref prog_zcv "zcv" program.

   If there is only one argument on the command line, @b zpr writes to stdout.
   A second argument, if present, is taken as the output file name.

   To find out the contents of a MeatAxe file, use the -s option. To generate
   output readble by GAP, use -G. Both options can be combined.
 */
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
