////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Test the arithmetic module.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"

#include <stdlib.h>
#include <string.h>

#include <test_table.c>

void test_Fail(const char *file, int line, const char *func, const char *msg, ...)
{
    va_list args;
    va_start(args, msg);
    fprintf(stderr, "%s:%d:: error: TEST FAILED: %s\n", file, line, func);
    vfprintf(stderr, msg, args);
    fprintf(stderr, "\n");
    exit(1);
}

void test_Assert(const char *file, int line, const char *func, int e, const char *estr)
{
   if (e) return;
   test_Fail(file, line, func, "assertion failed: %s", estr);
}

void test_EqInt(const char *file, int line, const char *func, int act, int exp,
	        const char *actstr, const char *expstr)
{
   if (act == exp) return;  
   test_Fail(file, line, func, "value of %s:\nactual:   %d\nexpected: %d (%s)\n",
	     actstr, act, exp, expstr);
}

static int Fields[] = {2,3,4,5,16,67,125,256,-1};
static int NextFieldIndex = 0;
FEL *FTab = NULL;
static MtxApplicationInfo_t AppInfo = {
   "mtxtest", "MeatAxe Library test program",
   "SYNTAX\n"
   "    mtxtest " MTX_COMMON_OPTIONS_SYNTAX " [-t <Field>]\n"
   "\n"
   "ARGUMENTS\n"
   "\n"
   "OPTIONS\n"
   "    -t <Field> .............. Print tables for GF(<Field>)\n"
   MTX_COMMON_OPTIONS_DESCRIPTION
};


/* --------------------------------------------------------------------------
   Error() - Print error message and exit
   -------------------------------------------------------------------------- */

void Error(char *msg, ...)
{
   va_list al;
   va_start(al,msg);
   fprintf(stderr,"\n*** ERROR:");
   vfprintf(stderr,msg,al);
   fprintf(stderr,"\n");
   va_end(al);
   exit(1);
}


/* --------------------------------------------------------------------------
   MakeFTab() - Create an array of all field elements.
   -------------------------------------------------------------------------- */

void MakeFTab()
{
   int i;

   if (FTab != NULL) {
      free(FTab);
   }
   FTab = NALLOC(FEL,FfOrder);
   for (i = 0; i < FfOrder; ++i) {
      FTab[i] = FfFromInt(i);
      if (!ISFEL(FTab[i])) {
         Error("FfFromInt(%d)=%d, illegal value",i,FTab[i]);
      }
   }
}


/* --------------------------------------------------------------------------
   NextField() - Set the next field
   -------------------------------------------------------------------------- */

int NextField()
{
   int f = Fields[NextFieldIndex];
   if (f > 0) {
      ++NextFieldIndex;
      FfSetField(f);
      MakeFTab();
   }
   return f;
}


/* --------------------------------------------------------------------------
   SelectField() - Set a specific field
   -------------------------------------------------------------------------- */

void SelectField(int f)
{
   FfSetField(f);
   MakeFTab();
}


/* --------------------------------------------------------------------------
   MkMat() - Create a matrix
   -------------------------------------------------------------------------- */

Matrix_t *MkMat(int nor, int noc, ...)
{
   int i;
   Matrix_t *m;
   va_list al;
   va_start(al,noc);

   m = MatAlloc(FfOrder,nor,noc);
   for (i = 0; i < nor; ++i) {
      PTR x = MatGetPtr(m,i);
      int k;
      for (k = 0; k < noc; ++k) {
         int a = va_arg(al,int);
         if (a >= 0) {
            FfInsert(x,k,FTab[a % FfChar]);
         } else {
            FfInsert(x,k,FfNeg(FTab[-a % FfChar]));
         }
      }
   }
   return m;
}


/* --------------------------------------------------------------------------
   Test function table
   -------------------------------------------------------------------------- */

#define CHECK_FUNCTION_TABLE

static int prtables(int field)
{
   int a, b;

   FfSetField(field);
   printf(" + ");
   for (a = 0; a < FfOrder; ++a) {
      printf("%3d", a);
   }
   printf("\n");
   for (a = 0; a < FfOrder; ++a) {
      printf("%3d",a);
      for (b = 0; b < FfOrder; ++b) {
         printf("%3d",FfToInt(FfAdd(FfFromInt(a),FfFromInt(b))));
      }
      printf("\n");
   }

   printf("\n * ");
   for (a = 0; a < FfOrder; ++a) {
      printf("%3d", a);
   }
   printf("\n");
   for (a = 0; a < FfOrder; ++a) {
      printf("%3d",a);
      for (b = 0; b < FfOrder; ++b) {
         printf("%3d",FfToInt(FfMul(FfFromInt(a),FfFromInt(b))));
      }
      printf("\n");
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static int TstErrorHandlerActive = 0;
static int TstErrorCode = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TstErrorHandler(const MtxErrorRecord_t *err)
{
   TstErrorCode = (err != NULL) ? 1 : 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void TstClearError()
{
   TstErrorCode = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void TstStartErrorChecking()
{
   if (TstErrorHandlerActive) {
       Error("TstStartErrorChecking(): already started");
   }
   if (MtxSetErrorHandler(TstErrorHandler)) {
       Error("TstStartErrorChecking(): other handler is set");
   }
   TstErrorHandlerActive = 1;
   TstClearError();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void TstStopErrorChecking()
{
   if (TstErrorHandlerActive) {
      TstClearError();
      MtxSetErrorHandler(NULL);
      TstErrorHandlerActive = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns and clears the error status.

int TstHasError()
{
   const int ec = TstErrorCode;
   TstErrorCode = 0;
   return ec != 0;
}



/* --------------------------------------------------------------------------
   main()
   -------------------------------------------------------------------------- */

int main(int argc, const char **argv)
{
   MtxApplication_t *app;
   int i;
   int field;

   if ((app = AppAlloc(&AppInfo,argc,argv)) == NULL) {
      return -1;
   }
   field = AppGetIntOption(app,"-t --print-tables",-1,2,256);
   if (field > 0) {
      prtables(field);
      exit(0);
   }

#if 0
   SelectField(2);
   {
      Matrix_t *mat = RndMat(2,1000,1000);
      Matrix_t *a = RndMat(2,1,1000);
      Matrix_t *b = RndMat(2,1,1000);
      PTR mm = mat->Data, aa = a->Data, bb = b->Data;
      for (i = 0; i < 50000; ++i) {
         FfMapRow(aa,mm,1000,bb);
      }
      return 0;
   }
#endif

#if 0
   SelectField(2);
   {
      Matrix_t *mat = RndMat(2,1000,1000);
      Matrix_t *a = RndMat(2,1,1000);
      Matrix_t *b = RndMat(2,1,1000);
      GreasedMatrix_t *gm = GrMatAlloc(mat,4);
      PTR aa = a->Data, bb = b->Data;
      for (i = 0; i < 50000; ++i) {
         GrMapRow(aa,gm,bb);
      }
      return 0;
   }
#endif

#if 0
   {
      Matrix_t *id = MatId(3,5);
      GreasedMatrix_t *gm = GrMatAlloc(id,2);
      MtxFile_t *f = MfCreate("greas.test",gm->Field,gm->NumVecs,gm->Noc);
      MfWriteRows(f,gm->PrecalcData,gm->NumVecs);
      MfClose(f);
   }
#endif

   printf("MeatAxe Version %s\n",MtxVersion);

   for( i = 0; i < sizeof(test_AllTests)/sizeof(test_AllTests[0]); ++i) {
       const test_Definition * const td = test_AllTests + i;
       printf("+ %s\n", td->name);
       fflush(stdout);
       test_AllTests[i].f();
       if (TstHasError()) {
          test_Fail(td->file,td->line,td->name,"Uncaught error");
       }
       TstStopErrorChecking();
       NextFieldIndex = 0;
   }
   printf("All tests passed\n");
   return 0;
}
