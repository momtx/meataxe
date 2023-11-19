////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Playground
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

/* ------------------------------------------------------------------
   Variables
   ------------------------------------------------------------------ */

static MtxApplicationInfo_t AppInfo = {
   "playground", "XXX",
   "SYNTAX\n"
   "    playground " MTX_COMMON_OPTIONS_SYNTAX "XXX\n"
   "\n"
   "ARGUMENTS\n"
   "    XXX ..................... ...\n"
   "\n"
   "OPTIONS\n"
   MTX_COMMON_OPTIONS_DESCRIPTION
   "\n"
   "FILES\n"
   "    XXX ..................... I XXX\n"
   "    XXX ..................... O XXX\n"
};

static MtxApplication_t* App = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////

//static uint32_t nextPrime(uint32_t q)
//{
//   if (q < 2) return 2;
//   if (q == 2) return 3;
//
//   if (q % 2 == 0) ++q; else q += 2;
//   while (1) {
//      for (uint32_t i = 3; ; i += 2) {
//         if ((q % i) == 0) break;
//         if (i * i >= q)
//            return q;
//      }
//      q += 2;
//   }
//}

//static uint32_t nextPrimePower(uint32_t q)
//{
//   if (q < 2) return 2;
//   if (q == 2) return 3;
//
//   while (1) {
//      ++q;
//      if (q % 2 == 0) {
//         uint32_t r;
//         for (r = q; r != 0 && r % 2 == 0; r /= 2);
//         if (r == 1) return q;
//      } else {
//         for (uint32_t i = 3; ; i += 2) {
//            if ((q % i) == 0) break;
//            if (i * i >= q)
//               return q;
//         }
//      }
//   }
//}


//static void test()
//{
//   printf("ZZZ=%d\n", MTX_ZZZ);
//   for (uint32_t q = 2; q <= 256; q = nextPrimePower(q))
//   {
//      ffSetField(q);
//      printf("q=%4u   gen=0x%02x (%3d)\n", (unsigned) q, (unsigned) ffGen, ffToInt(ffGen));
//   }
//}

static void test()
{
   char* tmp[25];
   const size_t N = sizeof(tmp) / sizeof(tmp[0]);
   for (int i = 0; i < N; ++i)
   {
      StrBuffer* sb = sbAlloc(10);
      sbPrintf(sb, "test %d/%d ___________________________________________________", i, (int)N);
      tmp[i] = sbToEphemeralString(sb);
   }
   for (int i = 0; i < N; ++i)
   {
      printf("check %d: %s\n", i, tmp[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
   App = appAlloc(&AppInfo, argc, argv);
   appGetArguments(App, 0, 0);
   test();
   return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
