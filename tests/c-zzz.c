////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Test the arithmetic module.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <test_table.c>


void tstPrintRows(const char *name, PTR x, int nor, int noc)
{
   printf("---\n%s (%dx%d):\n", name, nor, noc);
   for (int n = nor; n > 0; --n) {
      for (int col = 0; col < noc; ++col) {
         printf(" %5d",ffToInt(ffExtract(x,col)));
      }
      printf("\n");
      x = (PTR)((char*) x + ffRowSize(noc));
   }
}

static char *argv0;
static const char* tstCurrent = "";
static int tstFailCalled = 0;
static int tstMessageThreshold = 0;


static void tstVprintf(int level, const char *msg, va_list args)
{
   if (level <= tstMessageThreshold)
      vfprintf(stdout, msg, args);
}

MTX_PRINTF(2,3)
static void tstPrintf(int level, const char *msg, ...)
{
   va_list args;
   va_start(args, msg);
   tstVprintf(level, msg, args);
   va_end(args);
}


void tstFail(const char *file, int line, const char *func, const char *msg, ...)
{
    va_list args;
    va_start(args, msg);
    if (strcmp(tstCurrent, func) == 0) {
       tstPrintf(-1,"%s:%d:: error: TEST FAILED: %s\n", file, line, tstCurrent);
    } else {
       tstPrintf(-1,"%s:%d:: error: TEST FAILED: %s (%s)\n", file, line, tstCurrent, func);
    }
    tstVprintf(-1,msg, args);
    tstPrintf(-1,"\n");
    tstFailCalled = 1;
}

int tstAssert(const char *file, int line, const char *func, int e, const char *estr)
{
   if (e) return 0;
   tstFail(file, line, func, "assertion failed: %s", estr);
   return 1;
}

int tstAssertEqInt(const char *file, int line, const char *func, int act, int exp,
	        const char *actstr, const char *expstr)
{
   if (act == exp) return 0;  
   tstFail(file, line, func, "value of %s:\nactual:   %d\nexpected: %d (%s)\n",
	     actstr, act, exp, expstr);
   return 1;
}

static const int DEFAULT_FIELDS[] = {2,3,4,5,16,67,125,256,
#if MTX_ZZZ == 1
    59049, // = 3 ^ 10
#endif
    -1};

static const int* SelectedFields = DEFAULT_FIELDS;
static const int* CurrentField = NULL;
static int defaultField = 243;



FEL *FTab = NULL;
static MtxApplicationInfo_t AppInfo = {
   "mtxtest", "MeatAxe Library test program",
   "SYNTAX\n"
   "    mtxtest " MTX_COMMON_OPTIONS_SYNTAX " [-l] [-t <Field>] [-f <Field>] [<TestSpec>]\n"
   "\n"
   "ARGUMENTS\n"
   "    <TestSpec> .............. Test(s) to be run (shell-style pattern)\n"
   "\n"
   "OPTIONS\n"
   "    -t, --print-tables ...... Print tables for GF(<Field>)\n"
   "    -f, --field ............. Execute tests only for a single field.\n"
   "    -l, --list-tests ........ List all avaliable tests and exit\n"
   MTX_COMMON_OPTIONS_DESCRIPTION
};


/* --------------------------------------------------------------------------
   MakeFTab() - Create an array of all field elements.
   -------------------------------------------------------------------------- */

void MakeFTab()
{
   free(FTab);
   FTab = NALLOC(FEL, ffOrder);
   for (int i = 0; i < ffOrder; ++i) {
      FTab[i] = ffFromInt(i);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

void SelectField(int f)
{
   ffSetField(f);
   MakeFTab();
   RngReset();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int NextField()
{
   if (CurrentField == NULL)
      CurrentField = SelectedFields;
   else if (*CurrentField != -1)
      ++CurrentField;
   if (*CurrentField == -1)
      return 0;
   SelectField(*CurrentField);
   ffSetField(*CurrentField);
   MakeFTab();
   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void UseFixedField(int field)
{
   static int fixedField[2];
   fixedField[0] = field;
   fixedField[1] = -1;
   SelectedFields = fixedField;
   defaultField = field;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static uint32_t rng = 0;
void RngReset()
{
    rng = 0;
}

uint32_t RngNext()
{
    rng = rng * 69069U + 107U;
    return rng;
}

FEL RandomFieldElement()
{
    return FTab[RngNext() % ffOrder];
}

FEL RandomNonzeroFieldElement()
{
    while (1) {
	FEL x = FTab[RngNext() % ffOrder];
	if (x != FF_ZERO)
	    return x;
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

void ForEachField(const char *testName, void (*testFunction)())
{
   for (const int *optr = SelectedFields; *optr > 0; ++optr) {
      printf("+ %s - GF(%d)\n", testName, *optr);
      SelectField(*optr);
      RngReset();
      testFunction();
   }
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

   m = matAlloc(ffOrder,nor,noc);
   for (i = 0; i < nor; ++i) {
      PTR x = matGetPtr(m,i);
      int k;
      for (k = 0; k < noc; ++k) {
         int a = va_arg(al,int);
         if (a >= 0) {
            ffInsert(x,k,FTab[a % ffChar]);
         } else {
            ffInsert(x,k,ffNeg(FTab[-a % ffChar]));
         }
      }
   }
   return m;
}


/* --------------------------------------------------------------------------
   Test function table
   -------------------------------------------------------------------------- */

#define CHECK_FUNCTION_TABLE

static void printTables(int field)
{
   int a, b;

   int width = field <= 256 ? 3 : 6;

   ffSetField(field);
   printf("%*s + ", width - 3, "");
   for (a = 0; a < ffOrder; ++a) {
      printf("%*d", width, a);
   }
   printf("\n");
   for (a = 0; a < ffOrder; ++a) {
      printf("%*d", width,a);
      for (b = 0; b < ffOrder; ++b) {
         printf("%*d", width,ffToInt(ffAdd(ffFromInt(a),ffFromInt(b))));
      }
      printf("\n");
   }

   printf("\n");
   printf("%*s * ", width - 3, "");
   for (a = 0; a < ffOrder; ++a) {
      printf("%*d", width, a);
   }
   printf("\n");
   for (a = 0; a < ffOrder; ++a) {
      printf("%*d", width,a);
      for (b = 0; b < ffOrder; ++b) {
         printf("%*d",width, ffToInt(ffMul(ffFromInt(a),ffFromInt(b))));
      }
      printf("\n");
   }

   printf("\n");
   printf("Subfield elements\n");
   const size_t MAX_EMBED = sizeof(mtx_subfields) / sizeof(mtx_subfields[0]);

   for (size_t i = 0; i < MAX_EMBED && mtx_subfields[i] >= 2; ++i) {
      const int subfield = mtx_subfields[i];
      // Get all subfield elements (temporary swith the working field).
      FEL* subfieldElements = NALLOC(FEL, subfield);
      ffSetField(subfield);
      for (int k = 0; k < subfield; ++k)
         subfieldElements[k] = ffFromInt(k);

      // Print the embedding.
      ffSetField(field);
      printf("%5d   ", subfield);
      for (int k = 0; k < subfield; ++k) {
         printf(" %d", ffToInt(ffEmbed(subfieldElements[k],subfield)));
      }
      printf("\n");
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

struct TstAbortState tstAbortState = {0};


static void CatchAbortHandler(const struct MtxErrorInfo *err)
{
   if (tstAbortState.enabled) {
      longjmp(tstAbortState.jumpTarget, 112);
   }
   tstFail(__FILE__, __LINE__, __func__,
         "UNEXPECTED ABORT\nabort reason: %s\nCANNOT CONTINUE TESTS, EXITING",
         err->message);
   exit(2);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void tstPrepareCatchAbort(const char *file, int line, const char *func, const char *estr)
{
   tstAbortState.file = file;
   tstAbortState.line = line;
   tstAbortState.func = func;
   tstAbortState.expr = estr;
   tstAbortState.enabled = 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void tstMissingAbort()
{
   tstAbortState.enabled = 0;
   tstFail(tstAbortState.file, tstAbortState.line, tstAbortState.func,
         "Did not abort as expected\nexpr: %s", tstAbortState.expr);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int strmatch(const char *s, const char *pattern)
{
   while (*pattern) {
      if (*pattern == '?') {
         if (*s == 0) return 1;
         ++s;
      } else if (*pattern == '*') {
         while (1) {
            if (!strmatch(s,pattern+1)) return 0;
            if (*s == 0) return 1;
            ++s;
         }
      } else {
         if (*s != *pattern) return 1;
         ++s;
      }
      ++pattern;
   }
   return *s == 0 ? 0 : 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int testSelected(const char *name, int nsel, const char * const *sel)
{
   if (nsel == 0) return 1;	// no selection, run all tests
   for (int i = 0; i < nsel; ++i) {
      if (strmatch(name,sel[i]) == 0) return 1;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

typedef TstResult (*SimpleTestFunction)(void);
typedef TstResult (*FieldDependentTestFunction)(int q);

static int executeTest(const struct TstFoundTest* test, int field)
{
   mtxCleanupLibrary();
   unsetenv("MTXLIB");
   mtxInitLibrary(argv0);
   mtxRandomInit(52134);
   TstResult result = 0;
   tstFailCalled = 0;

   MtxSetErrorHandler(CatchAbortHandler);
   tstAbortState.enabled = 0;
   tstCurrent = test->name;

   if (field > 0) {
      tstPrintf(0,"+ %s - GF(%d)\n", test->name, field);
      SelectField(field);
      FieldDependentTestFunction tf = (FieldDependentTestFunction) test->f;
      result = tf(field);
   } else {
      tstPrintf(0,"+ %s\n", test->name);
      SelectField(defaultField);
      SimpleTestFunction tf = (SimpleTestFunction) test->f;
      result = tf();
   }
   if (result != 0 && !tstFailCalled)
      tstFail(__FILE__, __LINE__, __func__, "Test failed with no error message");
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void listTests(int nsel, const char* const* sel)
{
   for (struct TstFoundTest* t = foundTests; t->f != NULL; ++t) {
      if (!testSelected(t->name, nsel, sel)) continue;
      printf("%s%s\n", t->name, (t->flags & TST_FLAG_PER_FIELD) ? "(q)" : "");
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int cmpTests(const void* a, const void* b)
{
   const struct TstFoundTest* ta = (const struct TstFoundTest*)a;
   const struct TstFoundTest* tb = (const struct TstFoundTest*)b;
   return strcmp(ta->name, tb->name);
}

static void sortTests()
{
   qsort(foundTests, sizeof(foundTests) / sizeof(foundTests[0]) - 1,
         sizeof(foundTests[0]), cmpTests);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   argv0 = strdup(argv[0]);
   MtxApplication_t *app;

   if ((app = appAlloc(&AppInfo,argc,argv)) == NULL) {
      return -1;
   }
#if MTX_ZZZ == 1
   static const int MTX_MAX_Q = 65535;
#elif MTX_ZZZ == 0
   static const int MTX_MAX_Q = 256;
#else
#error MTX_ZZZ undefined
#endif

   int field = appGetIntOption(app,"-t --print-tables",MTX_NVAL, 2, MTX_MAX_Q);
   if (field != MTX_NVAL) {
      printTables(field);
      exit(0);
   }

   field = appGetIntOption(app,"-f --field",MTX_NVAL, 2, MTX_MAX_Q);
   if (field != MTX_NVAL) {
      UseFixedField(field);
   }
   const int listOnly = appGetOption(app,"-l --list-tests" );
   tstMessageThreshold = MtxMessageLevel;
   MtxMessageLevel = 0;

   const int nsel = appGetArguments(app,0,1000);
   const char * const * const sel = (const char* const*)app->argV;

   sortTests();

   if (listOnly) {
      listTests(nsel, sel);
      return 0;
   }

   tstPrintf(0,"MeatAxe Version %s\n",mtxVersion());

   int nAvailable = 0;
   int nSelected = 0;
   int nFailed = 0;

   // Execute field-dependent tests
   for (const int* currentField = SelectedFields; *currentField > 1; ++currentField) {
      for (const struct TstFoundTest* t = foundTests; t->f != NULL; ++t) {
         if ((t->flags & TST_FLAG_PER_FIELD) == 0) continue;
         ++nAvailable;
         if (!testSelected(t->name, nsel, sel)) continue;
         ++nSelected;
         if (executeTest(t, *currentField) != 0)
            ++nFailed;
      }
   }

   // Execute field-independent tests
   for (struct TstFoundTest* t = foundTests; t->f != NULL; ++t) {
      if ((t->flags & TST_FLAG_PER_FIELD) != 0) continue;
      ++nAvailable;
      if (!testSelected(t->name, nsel, sel)) continue;
      ++nSelected;
      if (executeTest(t, -1) != 0)
         ++nFailed;
   }

   tstPrintf(-2,"\nTest results: %d total, %d selected", nAvailable, nSelected);
   if (nFailed == 0)
      tstPrintf(-2, " -- no failures\n");
   else
      tstPrintf(-2, ", %d FAILED\n", nFailed);
   return nFailed > 0 ? 1 : 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
