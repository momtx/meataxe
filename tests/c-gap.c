////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for GAP formatting functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>


static int assertGapFormat(
      const struct TstSourceLocation* where, FEL f, const char *fExpr, const char *expected)
{
   const char *actual = gapFelToString(f);
   if (strcmp(actual, expected) == 0)
      return 0;
   tstFail(where, "Wrong GAP representation of 0x%04x (%s):\nactual:   \"%s\"\nexpected: \"%s\"\n",
         (unsigned) f, fExpr, actual, expected);
   return 1;

}

#define ASSERT_GAP_FORMAT(expr, fmt, ...) \
   do {\
      char tmp[100]; \
      snprintf(tmp, sizeof(tmp), fmt, __VA_ARGS__); \
      if (assertGapFormat(TST_HERE, expr, #expr, tmp) != 0) return 1; \
   } while (0)

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Gap_FfGenerator(int q)
{
   if (ffChar == ffOrder) {
      ASSERT_GAP_FORMAT(ffGen, "Z(%u)*1", ffOrder);
   } else {
      ASSERT_GAP_FORMAT(ffGen, "Z(%u)^1", ffOrder);
   }
   return 0;
}

TstResult Gap_FfZero(int q)
{
   ASSERT_GAP_FORMAT(FF_ZERO, "Z(%u)*0", ffOrder);
   return 0;
}

TstResult Gap_FfOne(int q)
{
   if (ffChar == ffOrder) {
      unsigned n = 1;
      for (FEL a = ffGen; a != FF_ONE; a = ffAdd(a, ffGen)) {
         ++n;
      }
      ASSERT_GAP_FORMAT(FF_ONE, "Z(%u)*%u", ffOrder, n);
   } else {
      ASSERT_GAP_FORMAT(FF_ONE, "Z(%u)^0", ffOrder);
   }
   return 0;
}

TstResult Gap_PrimeFieldElements(int q)
{
   // note: i = 0,1 are checked in separate tests
   for (int i = 2; i < ffChar; ++i) {
      const FEL a = ffFromInt(i);
      if (ffChar == ffOrder) {
         unsigned n = 1;
         for (FEL b = ffGen; b != a; b = ffAdd(b, ffGen)) {
            ++n;
         }
         ASSERT_GAP_FORMAT(a, "Z(%u)*%u", ffOrder, n);


         char buf[20];
         snprintf(buf, sizeof(buf), gapFelToString(a));
         char *star = strchr(buf, '*');
         ASSERT(star != NULL);
         *star = 0;
         const char *str1Expected = star + 1;
         const char *str2Expected = buf;
         ASSERT_EQ_STRING(gapFelToString1(a), str1Expected);
         ASSERT_EQ_STRING(gapFelToString2(), str2Expected);
      }
      else {
         unsigned n = 1;
         for (FEL b = ffGen; b != a; b = ffMul(b, ffGen)) {
            ++n;
         }
         ASSERT_GAP_FORMAT(a, "Z(%u)^%u", ffOrder, n);
      }
   }
   return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
