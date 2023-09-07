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

#define ISFEL(f) ((f) == 0xffff || (unsigned short)(f) < (unsigned short)ffOrder)

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

int tstAssert(const char *file, int line, const char *func, int e, const char *estr);
#define ASSERT(e) \
   do { if (tstAssert(__FILE__,__LINE__,__func__,e,#e) != 0) return 1; } while (0)

int tstAssertEqInt(const char *file, int line, const char *func, int act, int exp,
	        const char *actstr, const char *expstr);
#define ASSERT_EQ_INT(act,exp) \
   do { if (tstAssertEqInt(__FILE__,__LINE__,__func__,act,exp,#act,#exp) != 0) return 1; } while(0)

struct TstAbortState {
   jmp_buf jumpTarget;
   int enabled;
   const char* file;
   int line;
   const char* func;
   const char* expr;
};
extern struct TstAbortState tstAbortState;
void tstPrepareCatchAbort(const char *file, int line, const char *func, const char *estr);
void tstMissingAbort();
#define ASSERT_ABORT(e) \
   do {\
      tstPrepareCatchAbort(__FILE__, __LINE__, __func__, #e); \
      if (setjmp(tstAbortState.jumpTarget) == 0) { \
         e; \
         tstMissingAbort(); \
      } \
      tstAbortState.enabled = 0; \
   } while (0)


MTX_PRINTF_ATTRIBUTE(4,5)
void tstFail(const char *file, int line, const char *func, const char *msg, ...);
#define TST_FAIL(msg, ...) \
   do { tstFail(__FILE__, __LINE__, __func__, msg, __VA_ARGS__); return 1; } while(0)

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
