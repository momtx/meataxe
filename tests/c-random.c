////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for the pseudorandom number generator
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <string.h>


////////////////////////////////////////////////////////////////////////////////////////////////////

static int Test11(unsigned seed, const long *table)
{
   int i;

   mtxRandomInit(seed);
   for (i = 0; i < 10; ++i, ++table) {
      int k;
      long val = mtxRandom() & 0x7FFFFFFF;
      ASSERT(val == *table);
      for (k = 0; k < 61; ++k) {
         mtxRandom();
      }
   }
   return 0;
}

TstResult RandomNumberGenerator1()
{
   static long Table0[10] = {
      826837439, 1433481918, 1807203728, 1251143873,
      498964889,886423565,167672701, 1728315981,248403305,1037767977
   };
   static long Table1[10] = {
      269167349, 1677366103,1597714250, 970290675,
      436236141, 2108708678, 89648197, 1313827126, 514978688, 628812726
   };
   int result = 0;
   result |= Test11(0,Table0);
   result |= Test11(1,Table1);
   result |= Test11(0,Table0);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Check if the random numbers are sufficiently equally distributed

TstResult RandomNumberGenerator2()
{
   int count[100];
   int n;
   for (n = 10; n < 100; ++n) {

       const int loop = 1550;
       const int min = loop - loop / 10;
       const int max = loop + loop / 10;

       memset(count, 0, sizeof(count));
       for (int i = 0; i < loop * n; ++i) {
	   ++count[mtxRandomInt(n)];
       }

       for (int i = 0; i < n; ++i) {
	   ASSERT(count[i] >= min && count[i] <= max);
       }

   }
   return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
