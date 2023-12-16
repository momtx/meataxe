////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for the word generator
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////


static MatRep_t* makeRep(int field, int ngen, uint32_t dim)
{
   MatRep_t* mr = mrAlloc(0, NULL, 0);
   for (uint32_t i = 0; i < ngen; ++i) {
      Matrix_t *g = matAlloc(field, dim, dim);
      for (uint32_t k = 0; k < dim; ++k) {
         ffInsert(matGetPtr(g, k), (i + k) % dim, FF_ONE);
         ffInsert(matGetPtr(g, k), (i * k) % dim, FF_ONE);
      }
      mrAddGenerator(mr, g, 0);
   }
   return mr;
}

static WgData_t* makeWgen(int field, int ngen, int dim)
{
    return wgAlloc(makeRep(field, ngen, dim));
}

static void destroy(WgData_t* wg)
{
   MatRep_t* rep = (MatRep_t*)wg->Rep;
   wgFree(wg);
   mrFree(rep);
}

TstResult WordGenerator_RejectsWordNumberZero(int q)
{
   WgData_t* wg = makeWgen(q, 2, 1);
   ASSERT_ABORT(wgMakeWord(wg, 0));
   ASSERT_ABORT(wgSymbolicName(wg, 0));
   destroy(wg);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int WordGenerator_SymbolicName_(WgData_t* wg)
{
   ASSERT_EQ_STRING(wgSymbolicName(wg, 1), "a+b+ab");
   ASSERT_EQ_STRING(wgSymbolicName(wg, 2), "a+b+ab+ab2");
   ASSERT_EQ_STRING(wgSymbolicName(wg, 3), "a+ba+b2+bab+bab2");
   ASSERT_EQ_STRING(wgSymbolicName(wg, 4), "a+b+ba+b2+bab+bab2");
   ASSERT_EQ_STRING(wgSymbolicName(wg, 5), "a+b+ab+ba+b2+bab+bab2");
   ASSERT_EQ_STRING(wgSymbolicName(wg, 6), "a+ba+b2+ab2+bab+bab2");
   ASSERT_EQ_STRING(wgSymbolicName(wg, 7), "ab2+bab+bab2");
   ASSERT_EQ_STRING(wgSymbolicName(wg, 8), "a+b");
   ASSERT_EQ_STRING(wgSymbolicName(wg, 9), "a+ab");
   ASSERT_EQ_STRING(wgSymbolicName(wg,10), "b+ab");
   ASSERT_EQ_STRING(wgSymbolicName(wg,11), "a+ba");
   ASSERT_EQ_STRING(wgSymbolicName(wg,12), "b+ba");
   ASSERT_EQ_STRING(wgSymbolicName(wg,13), "a+b+ba");
   ASSERT_EQ_STRING(wgSymbolicName(wg,14), "ab+ba");
   ASSERT_EQ_STRING(wgSymbolicName(wg,15), "a+ab+ba");
   ASSERT_EQ_STRING(wgSymbolicName(wg,16), "b+ab+ba");
   ASSERT_EQ_STRING(wgSymbolicName(wg,17), "a+b+ab+ba");
   ASSERT_EQ_STRING(wgSymbolicName(wg,100), "b+ab+ab2+bab");
   ASSERT_EQ_STRING(wgSymbolicName(wg,1000), "babab+ba3b+a3ba");
   ASSERT_EQ_STRING(wgSymbolicName(wg,10000), "babab+ab3+aba2b+b2ab+a2ba+a2bab");
   return 0;
}

TstResult WordGenerator_SymbolicName(int q)
{
   WgData_t* wg = makeWgen(q, 2, 1);
   int result = WordGenerator_SymbolicName_(wg);
   destroy(wg);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int WordGenerator_SymbolicName3Gen_(WgData_t* wg)
{
   ASSERT_EQ_STRING(wgSymbolicName(wg, 1), "a+b+ca");
   ASSERT_EQ_STRING(wgSymbolicName(wg, 2), "a+b+ca+acb");
   ASSERT_EQ_STRING(wgSymbolicName(wg, 3), "a+cb+ba+c2b+b2ac");
   ASSERT_EQ_STRING(wgSymbolicName(wg, 4), "a+b+cb+ba+c2b+b2ac");
   ASSERT_EQ_STRING(wgSymbolicName(wg, 5), "a+b+ca+cb+ba+c2b+b2ac");
   ASSERT_EQ_STRING(wgSymbolicName(wg, 6), "a+cb+ba+acb+c2b+b2ac");
   ASSERT_EQ_STRING(wgSymbolicName(wg, 7), "acb+c2b+b2ac");
   ASSERT_EQ_STRING(wgSymbolicName(wg, 8), "a+b");
   ASSERT_EQ_STRING(wgSymbolicName(wg, 9), "a+ca");
   ASSERT_EQ_STRING(wgSymbolicName(wg,10), "b+ca");
   ASSERT_EQ_STRING(wgSymbolicName(wg,11), "a+cb");
   ASSERT_EQ_STRING(wgSymbolicName(wg,12), "b+cb");
   ASSERT_EQ_STRING(wgSymbolicName(wg,13), "a+b+cb");
   ASSERT_EQ_STRING(wgSymbolicName(wg,14), "ca+cb");
   ASSERT_EQ_STRING(wgSymbolicName(wg,15), "a+ca+cb");
   ASSERT_EQ_STRING(wgSymbolicName(wg,16), "b+ca+cb");
   ASSERT_EQ_STRING(wgSymbolicName(wg,17), "a+b+ca+cb");
   ASSERT_EQ_STRING(wgSymbolicName(wg,100), "b+ca+acb+c2b");
   ASSERT_EQ_STRING(wgSymbolicName(wg,1000), "ac2ac+ac2ab+ca4");
   ASSERT_EQ_STRING(wgSymbolicName(wg,10000), "ba2bc+ba2b+ba4+bc2a+c4+c2ac2");
   return 0;
}

TstResult WordGenerator_SymbolicName3Gen(int q)
{
   WgData_t* wg = makeWgen(q, 3, 1);
   int result = WordGenerator_SymbolicName3Gen_(wg);
   destroy(wg);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int assertFingerprint(const struct TstSourceLocation* where, WgData_t* wg,
      uint32_t n0, uint32_t n1, uint32_t n2, uint32_t n3, uint32_t n4, uint32_t n5)
{
   uint32_t fp[6];
   wgMakeFingerPrint(wg, fp);
   if (fp[0] != n0 || fp[1] != n1 || fp[2] != n2 || fp[3] != n3 || fp[4] != n4 || fp[5] != n5) {
      tstFail(where, "wrong fingerprint:\n"
            "actual:   %lu %lu %lu %lu %lu %lu\n"
            "expected: %lu %lu %lu %lu %lu %lu\n",
            (unsigned long) fp[0], (unsigned long) fp[1], (unsigned long) fp[2],
            (unsigned long) fp[3], (unsigned long) fp[4], (unsigned long) fp[5],
            (unsigned long) n0, (unsigned long) n1, (unsigned long) n2,
            (unsigned long) n3, (unsigned long) n4, (unsigned long) n5);
      return 1;
   }
   return 0;
}

#define ASSERT_FINGERPRINT(where, wg, n1, n2, n3, n4, n5, n6) \
   do { if (assertFingerprint(where, wg, n1, n2, n3, n4, n5, n6) != 0) return 1; } while (0)

static int WordGenerator_Fingerprint_(
      const struct TstSourceLocation* where, int field, int ngen,
      uint32_t n0, uint32_t n1, uint32_t n2, uint32_t n3, uint32_t n4, uint32_t n5)
{
   ffSetField(field);
   WgData_t* wg = makeWgen(field, ngen, 19);
   ASSERT_FINGERPRINT(where, wg, n0,n1,n2,n3,n4,n5);
   destroy(wg);
   return 0;
}

TstResult WordGenerator_Fingerprint()
{
   int result = 0;
   result |= WordGenerator_Fingerprint_(TST_HERE, 2, 2, 1,0,1,0,0,0);
   result |= WordGenerator_Fingerprint_(TST_HERE, 3, 2, 0,1,0,1,1,0);
   result |= WordGenerator_Fingerprint_(TST_HERE,64, 2, 1,0,1,0,0,0);
   
   result |= WordGenerator_Fingerprint_(TST_HERE, 2, 3, 1,0,0,1,0,1);
   result |= WordGenerator_Fingerprint_(TST_HERE, 3, 3, 1,0,0,0,0,0);
   result |= WordGenerator_Fingerprint_(TST_HERE,64, 3, 1,0,0,1,0,1);
   return result;
}


// vim:fileencoding=utf8:sw=3:ts=8:et:cin
