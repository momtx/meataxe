////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for bit strings
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

//#include "c-bitstring.h"
#include "check.h"

#include <string.h>

MTX_DEFINE_FILE_INFO

static int ErrorFlag = 0;

static void MyErrorHandler(const MtxErrorRecord_t *err)
{
   err = NULL;
   ErrorFlag = 1;
}


#ifdef DEBUG

static int CheckError()
{
   int i = ErrorFlag;
   ErrorFlag = 0;
   return i;
}


#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

#define NMAT 5

test_F TestBsAlloc()
{
   static int bssize[NMAT] = { 0,1,10,100,1000 };
   BitString_t  *m[NMAT];
   MtxErrorHandler_t *old_err_handler;
   int i;

   for (i = 0; i < NMAT; ++i) {
      m[i] = BsAlloc(bssize[i]);
   }
   for (i = 0; i < NMAT; ++i) {
      int k;
      BsIsValid(m[i]);
      MTX_VERIFY(m[i]->Size == bssize[i]);
      for (k = 0; k < bssize[i]; ++k) {
         ASSERT(!BsTest(m[i],k));
      }
   }
   for (i = 0; i < NMAT; ++i) {
      ASSERT_EQ_INT(BsFree(m[i]), 0);
   }
   old_err_handler = MtxSetErrorHandler(MyErrorHandler);
#ifdef DEBUG
   for (i = 0; i < NMAT; ++i) {
      if (BsIsValid(m[i]) || !CheckError()) { Error("BsIsValid() failed"); }}
#endif
   MtxSetErrorHandler(old_err_handler);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestSetClear(int size, BitString_t *a)
{
   int i;

   for (i = 0; i < size; ++i) {
      int k;
      for (k = 0; k < i; ++k) {
         ASSERT(BsTest(a,k));
      }
      for (; k < size; ++k) {
         ASSERT(!BsTest(a,k));
      }
      BsSet(a,i);
   }
   for (i = 0; i < size; ++i) {
      int k;
      for (k = 0; k < i; ++k) {
         ASSERT(!BsTest(a,k));
      }
      for (; k < size; ++k) {
         ASSERT(BsTest(a,k));
      }
      BsClear(a,i);
   }

   for (i = 0; i < size; ++i) {
      BsSet(a,i);
   }
   ASSERT_EQ_INT(BsClearAll(a), 0);
   for (i = 0; i < size; ++i) {
      ASSERT(!BsTest(a,i));
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

test_F BsBasicOperations()
{
   const int size = 50;
   BitString_t *a;
   a = BsAlloc(size);
   TestSetClear(size,a);
   BsFree(a);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestCompare1(int size, BitString_t *a, BitString_t *b)
{
   int i;
   for (i = 0; i < size; ++i) {
      if (BsCompare(a,b) != 0) { Error("BsCompare(x,x) != 0"); }
      BsSet(a,i);
      if (BsCompare(a,b) == 0) {
         Error("BsCompare() did not detect difference in bit %d",i);
      }
      BsSet(b,i);
   }
   for (i = 0; i < size; ++i) {
      if (BsCompare(a,b) != 0) { Error("BsCompare(x,x) != 0"); }
      BsClear(a,i);
      if (BsCompare(a,b) == 0) {
         Error("BsCompare() did not detect difference in bit %d",i);
      }
      BsClear(b,i);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

void TestBsCompare(unsigned flags)
{
   const int size = 50;
   BitString_t *a, *b;
   a = BsAlloc(size);
   b = BsAlloc(size);
   TestCompare1(size,a,b);
   BsFree(b);
   BsFree(a);
   flags = 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestCopy1(int size, BitString_t *a, BitString_t *b)
{
   for (int i = 0; i < size; i += 5) {
      BsSet(a,i);
   }
   BsCopy(b,a);
   ASSERT(!BsCompare(a,b));
   BitString_t *c = BsDup(a);
   ASSERT(!BsCompare(a,c));
   BsFree(c);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

test_F BsCopyAndDup()
{
   const int size = 49;
   BitString_t *a, *b;
   a = BsAlloc(size);
   b = BsAlloc(size);
   TestCopy1(size,a,b);
   BsFree(b);
   BsFree(a);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static BitString_t *RndBs(int size)
{
   int i;
   BitString_t *bs = BsAlloc(size);
   for (i = 0; i < size; ++i) {
      if (MtxRandomInt(2) != 0) {
         BsSet(bs,i);
      }
   }
   return bs;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void CheckIo1(BitString_t **bs, int n)
{
   const char file_name[] = "check.1";
   FILE *f;
   int i;

   /* Write bit strings
      ----------------- */
   f = SysFopen(file_name,FM_CREATE);
   for (i = 0; i < n; ++i) {
      BsWrite(bs[i],f);
   }
   fclose(f);

   /* Read bit strings
      ---------------- */
   f = SysFopen(file_name,FM_READ);
   for (i = 0; i < n; ++i) {
      BitString_t *a = BsRead(f);
      if (BsCompare(a,bs[i]) != 0) {
         Error("Read error");
      }
   }
   fclose(f);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

test_F BsFileIo()
{
   BitString_t *a[10];
   int i;

   MtxRandomInit(1235);
   a[0] = BsAlloc(0);
   for (i = 1; i < 10; ++i) {
      a[i] = RndBs(MtxRandomInt(100));
   }
   CheckIo1(a,10);
   for (i = 1; i < 10; ++i) {
      BsFree(a[i]);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestAndOr1(int size)
{
   BitString_t *a = BsAlloc(size);
   BitString_t *b = BsAlloc(size);
   BitString_t *bsor, *bsand, *bsminus;
   int i;

   for (i = 0; i < size; ++i) {
      if (i % 3 == 0) { BsSet(a,i); }
      if (i % 5 == 1) { BsSet(b,i); }
   }

   bsor = BsDup(a);
   BsOr(bsor,b);
   bsand = BsDup(a);
   BsAnd(bsand,b);
   bsminus = BsDup(a);
   BsMinus(bsminus,b);

   for (i = 0; i < size; ++i) {
      ASSERT_EQ_INT(BsTest(bsor,i), i % 3 == 0 || i % 5 == 1);
      ASSERT_EQ_INT(BsTest(bsand,i), i % 3 == 0 && i % 5 == 1);
      ASSERT_EQ_INT(BsTest(bsminus,i), i % 3 == 0 && i % 5 != 1);
   }

   BsFree(bsor);
   BsFree(bsand);
   BsFree(bsminus);
   BsFree(b);
   BsFree(a);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

test_F BsAndOrMinusRandomTest()
{
   for (int i = 0; i < 150; ++i) {
      TestAndOr1(i);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void CheckCount1(int size, BitString_t *a, BitString_t *b)
{
   int i;
   int count;

   for (count = 0, i = 0; i < size; ++i) {
      if (BsTest(a,i) && BsTest(b,i)) {
         ++count;
      }
   }
   ASSERT_EQ_INT(BsIntersectionCount(a,b), count);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

test_F BsIntersectionCountRandomTest()
{
   int i;
   MtxRandomInit(42);
   for (i = 0; i < 100; ++i) {
      int size = MtxRandomInt(100);
      BitString_t *a = RndBs(size), *b = RndBs(size);
      CheckCount1(size,a,b);
      BsFree(a);
      BsFree(b);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

test_F BsIsSubRandomTest()
{
   int i;
   BitString_t *a, *b;

   MtxRandomInit(42);
   for (i = 0; i < 150; ++i) {
      const int size = i;
      int k;
      a = RndBs(size);
      b = BsAlloc(size);
      ASSERT(BsIsSub(b,a));
      BsCopy(b,a);
      ASSERT(BsIsSub(b,a));
      for (k = 0; k < size; ++k) {
         if (!BsTest(a,k)) {
            BsSet(b,k);
            ASSERT(!BsIsSub(b,a));
            BsClear(b,k);
         }
      }
      for (k = 0; k < size; ++k) {
         BsClear(b,k);
         ASSERT(BsIsSub(b,a));
      }
      BsFree(a);
      BsFree(b);
   }
}
