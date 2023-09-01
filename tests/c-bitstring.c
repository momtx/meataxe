////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for bit strings
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "testing.h"

#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_AllocFree()
{
#define N 5
   static int bssize[N] = { 0,1,10,100,1000 };
   BitString_t *m[N];
   int i;

   for (i = 0; i < N; ++i) {
      m[i] = bsAlloc(bssize[i]);
   }
   for (i = 0; i < N; ++i) {
      int k;
      ASSERT(bsIsValid(m[i]));
      ASSERT(m[i]->Size == bssize[i]);
      for (k = 0; k < bssize[i]; ++k) {
         ASSERT(!bsTest(m[i],k));
      }
   }
   for (i = 0; i < N; ++i) {
      ASSERT_EQ_INT(bsFree(m[i]), 0);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_AbortsOnDoubleFree()
{
   BitString_t *bs = bsAlloc(100);
   ASSERT(bsIsValid(bs));
   bsFree(bs);
   ASSERT(!bsIsValid(bs));
   ASSERT_ABORT(bsFree(bs));
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestSetClear(int size, BitString_t *a)
{
   int i;

   for (i = 0; i < size; ++i) {
      int k;
      for (k = 0; k < i; ++k) {
         ASSERT(bsTest(a,k));
      }
      for (; k < size; ++k) {
         ASSERT(!bsTest(a,k));
      }
      bsSet(a,i);
   }
   for (i = 0; i < size; ++i) {
      int k;
      for (k = 0; k < i; ++k) {
         ASSERT(!bsTest(a,k));
      }
      for (; k < size; ++k) {
         ASSERT(bsTest(a,k));
      }
      bsClear(a,i);
   }

   for (i = 0; i < size; ++i) {
      bsSet(a,i);
   }
   ASSERT_EQ_INT(bsClearAll(a), 0);
   for (i = 0; i < size; ++i) {
      ASSERT(!bsTest(a,i));
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_SetClear()
{
   const int size = 50;
   BitString_t *a;
   a = bsAlloc(size);
   int result = TestSetClear(size,a);
   bsFree(a);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestCompare1(int size, BitString_t *a, BitString_t *b)
{
   int i;
   for (i = 0; i < size; ++i) {
      ASSERT_EQ_INT(bsCompare(a,b),0);
      bsSet(a,i);
      ASSERT(bsCompare(a,b) != 0);
      bsSet(b,i);
   }
   for (i = 0; i < size; ++i) {
      ASSERT_EQ_INT(bsCompare(a,b),0);
      bsClear(a,i);
      ASSERT(bsCompare(a,b) != 0);
      bsClear(b,i);
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Compare()
{
   const int size = 50;
   BitString_t *a, *b;
   a = bsAlloc(size);
   b = bsAlloc(size);
   int result = TestCompare1(size,a,b);
   bsFree(b);
   bsFree(a);
   return result;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestCopy1(int size, BitString_t *a, BitString_t *b)
{
   for (int i = 0; i < size; i += 5) {
      bsSet(a,i);
   }
   bsCopy(b,a);
   ASSERT(!bsCompare(a,b));
   BitString_t *c = bsDup(a);
   ASSERT(!bsCompare(a,c));
   bsFree(c);
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Copy()
{
   const int size = 49;
   BitString_t *a, *b;
   a = bsAlloc(size);
   b = bsAlloc(size);
   int result = TestCopy1(size,a,b);
   bsFree(b);
   bsFree(a);
   return result;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static BitString_t *RndBs(int size)
{
   int i;
   BitString_t *bs = bsAlloc(size);
   for (i = 0; i < size; ++i) {
      if (mtxRandomInt(2) != 0) {
         bsSet(bs,i);
      }
   }
   return bs;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static int CheckIo1(FILE* f, BitString_t **bs, int n)
{
   // Write bit strings
   for (int i = 0; i < n; ++i) {
      bsWrite(bs[i],f);
   }

   // Read bit strings
   rewind(f);
   for (int i = 0; i < n; ++i) {
      BitString_t *a = bsRead(f);
      ASSERT_EQ_INT(bsCompare(a,bs[i]), 0);
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitStringFileIo()
{
   BitString_t *a[10];
   int i;

   a[0] = bsAlloc(0);
   for (i = 1; i < 10; ++i) {
      a[i] = RndBs(mtxRandomInt(100));
   }
   const char file_name[] = "check.1";
   FILE* f = fopen(file_name,"w+");
   int result = CheckIo1(f, a, 10);
   fclose(f);
   for (i = 1; i < 10; ++i) {
      bsFree(a[i]);
   }
   return result;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestAndOr1(int size, BitString_t* a, BitString_t* b)
{
   BitString_t *bsor, *bsand, *bsminus;
   int i;

   for (i = 0; i < size; ++i) {
      if (i % 3 == 0) { bsSet(a,i); }
      if (i % 5 == 1) { bsSet(b,i); }
   }

   bsor = bsDup(a);
   bsOr(bsor,b);
   bsand = bsDup(a);
   bsAnd(bsand,b);
   bsminus = bsDup(a);
   bsMinus(bsminus,b);

   for (i = 0; i < size; ++i) {
      ASSERT_EQ_INT(bsTest(bsor,i), i % 3 == 0 || i % 5 == 1);
      ASSERT_EQ_INT(bsTest(bsand,i), i % 3 == 0 && i % 5 == 1);
      ASSERT_EQ_INT(bsTest(bsminus,i), i % 3 == 0 && i % 5 != 1);
   }

   bsFree(bsor);
   bsFree(bsand);
   bsFree(bsminus);
   return 0;
}

TstResult BitStringAndOrMinus()
{
   int result = 0;
   for (int size = 0; size < 150; ++size) {
      BitString_t *a = bsAlloc(size);
      BitString_t *b = bsAlloc(size);
      result |= TestAndOr1(size, a, b);
      bsFree(b);
      bsFree(a);
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int CheckCount1(int size, BitString_t *a, BitString_t *b)
{
   int i;
   int count;

   for (count = 0, i = 0; i < size; ++i) {
      if (bsTest(a,i) && bsTest(b,i)) {
         ++count;
      }
   }
   ASSERT_EQ_INT(bsIntersectionCount(a,b), count);
   return 0;
}

TstResult BitString_IntersectionCount()
{
   int result = 0;
   for (int i = 0; result == 0 && i < 100; ++i) {
      int size = mtxRandomInt(100);
      BitString_t *a = RndBs(size), *b = RndBs(size);
      result |= CheckCount1(size,a,b);
      bsFree(a);
      bsFree(b);
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int testIsSubset(BitString_t *a, BitString_t *b, int size)
{
   ASSERT(bsIsSub(b,a));
   bsCopy(b,a);
   ASSERT(bsIsSub(b,a));
   for (int k = 0; k < size; ++k) {
      if (!bsTest(a,k)) {
         bsSet(b,k);
         ASSERT(!bsIsSub(b,a));
         bsClear(b,k);
      }
   }
   for (int k = 0; k < size; ++k) {
      bsClear(b,k);
      ASSERT(bsIsSub(b,a));
   }
   return 0;
}

TstResult BitString_IsSubset()
{
   int i;
   BitString_t *a, *b;

   int result = 0;
   for (int size = 0; result == 0 && i < 150; ++i) {
      a = RndBs(size);
      b = bsAlloc(size);
      result = testIsSubset(a,b,size);
      bsFree(a);
      bsFree(b);
   }
   return result;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
