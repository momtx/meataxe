////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests framework
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _CHECK_H_
#define _CHECK_H_

#include <meataxe.h>

#include <stdarg.h>
#include <setjmp.h>

#if MTX_ZZZ == 0

#define ISFEL(f) ((unsigned int)(f) < (unsigned int)ffOrder)

#elif MTX_ZZZ == 1

#define ISFEL(f) ((f) == 0xffff || (uint16_t)(f) < (uint16_t)(ffOrder - 1))

#elif 
#error MTX_ZZZ undefined
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////////////////////////////////

void tstPrintRows(const char *name, PTR x, int nor, int noc);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Test framework
////////////////////////////////////////////////////////////////////////////////////////////////////

struct TstSourceLocation {
   const char* file;    ///< The source file name.
   int line;            ///< The line number.
   const char* func;    ///< The function name.
};

#define TST_HERE (&(const struct TstSourceLocation){__FILE__, __LINE__, __func__}) 

int tstAssert(const struct TstSourceLocation* where, int e, const char *estr);
#define ASSERT(e) \
   do { if (tstAssert(TST_HERE, e, #e) != 0) return 1; } while (0)

int tstAssertEqInt(const struct TstSourceLocation* where, int act, int exp,
	        const char *actstr, const char *expstr);
#define ASSERT_EQ_INT(act,exp) \
   do { if (tstAssertEqInt(TST_HERE,act,exp,#act,#exp) != 0) return 1; } while(0)

int tstAssertEqString(const struct TstSourceLocation* where, const char* act, const char* exp,
	        const char *actstr, const char *expstr);
#define ASSERT_EQ_STRING(act,exp) \
   do { if (tstAssertEqString(TST_HERE,act,exp,#act,#exp) != 0) return 1; } while(0)


struct TstAbortState {
   jmp_buf jumpTarget;
   int enabled;
   struct TstSourceLocation where;
   const char* expr;
};
extern struct TstAbortState tstAbortState;
void tstPrepareCatchAbort(const struct TstSourceLocation* where, const char *estr);
void tstMissingAbort();
#define ASSERT_ABORT(e) \
   do {\
      tstPrepareCatchAbort(TST_HERE, #e); \
      if (setjmp(tstAbortState.jumpTarget) == 0) { \
         e; \
         tstMissingAbort(); \
      } \
      tstAbortState.enabled = 0; \
   } while (0)


MTX_PRINTF(2,3)
void tstFail(const struct TstSourceLocation *where, const char *msg, ...);
#define TST_FAIL(msg, ...) \
   do { tstFail(TST_HERE, msg, __VA_ARGS__); return 1; } while(0)

extern FEL *FTab;
int NextField();
void SelectField(int f);
void MakeFTab();
void ForEachField(const char *testName, void (*testFunction)());

void RngReset();
uint32_t RngNext();
FEL RandomFieldElement();
FEL RandomNonzeroFieldElement();

extern Matrix_t *MkMat(int nor, int noc, ...);

Perm_t *RndPerm(int degree);
Matrix_t *RndMat(int fl, int nor, int noc);
Poly_t *RndPol(int fl, int mindeg, int maxdeg);

typedef int TstResult;

#define TST_FLAG_PER_FIELD 0x0001

struct TstFoundTest {
   void *f;
   unsigned flags;
   const char *name;
};



#endif

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
