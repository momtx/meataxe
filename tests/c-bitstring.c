////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for bit strings
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "testing.h"

#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_NullIsInvalid()
{
   ASSERT(!bsIsValid(NULL));
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

TstResult BitString_Fixed_AllocFree()
{
   static int SIZES[] = { 0,1,10,100, 1000,10000, 100000, -1 };
   for (int* size = SIZES; *size >= 0; ++size) {
      BitString_t* bs = bsAlloc(*size);
      ASSERT(bsIsValid(bs));
      ASSERT_EQ_INT(bs->size, *size);
      for (size_t i = 0; i < *size; ++i) {
         ASSERT(!bsTest(bs, i));
      }
      bsFree(bs);
      ASSERT(!bsIsValid(bs));
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Variable_AllocFree()
{
   BitString_t* bs = bsAllocEmpty();
   ASSERT(bsIsValid(bs));
   ASSERT_EQ_INT(bs->capacity, 0);
   bsFree(bs);
   ASSERT(!bsIsValid(bs));
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Variable_ReadingBeyondCapacity()
{
   BitString_t *bs = bsAllocEmpty();
   for (size_t i = 0; i < 10000; ++i)
      ASSERT_EQ_INT(bsTest(bs, i), 0);
   ASSERT_EQ_INT(bs->capacity, 0); // not extended by read
   bsFree(bs);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Variable_WritingBeyondCapacity()
{
   const size_t SIZE = 100;
   BitString_t* a1 = bsAllocEmpty();
   for (size_t i = 0; i < SIZE; i += 3)
      bsSet(a1, i);

   for (size_t i = 0; i < SIZE; ++i) {
      ASSERT_EQ_INT(bsTest(a1, i), (i % 3 == 0) ? 1 : 0);
   }
   for (size_t i = SIZE; i < 2*SIZE; ++i) {
      ASSERT_EQ_INT(bsTest(a1, i), 0);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Fixed_AbortsOnReadingBeyondEnd()
{
   static int SIZES[] = { 0,1,10,100, 1000,10000, 100000, -1 };
   for (int* size = SIZES; *size >= 0; ++size) {
      BitString_t* bs = bsAlloc(*size);
      ASSERT_ABORT(bsTest(bs, *size));
      bsFree(bs);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Fixed_AbortsOnWritingBeyondEnd()
{
   static int SIZES[] = { 0,1,10,100, 1000,10000, 100000, -1 };
   for (int* size = SIZES; *size >= 0; ++size) {
      BitString_t* bs = bsAlloc(*size);
      ASSERT_ABORT(bsSet(bs, *size));
      ASSERT_ABORT(bsClear(bs, *size));
      bsFree(bs);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Fixed_ClearAll()
{
   static size_t SIZE = 100;
   BitString_t* bs = bsAlloc(SIZE);
   for (size_t i = 0; i < SIZE; ++i)
      bsSet(bs, i);
   bsClearAll(bs);
   ASSERT_EQ_INT(bs->size, SIZE);
   for (size_t i = 0; i < SIZE; ++i)
      ASSERT_EQ_INT(bsTest(bs, i), 0);
   bsFree(bs);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Variable_ClearAll()
{
   static size_t SIZE = 100;
   BitString_t* bs = bsAllocEmpty();
   for (size_t i = 0; i < SIZE; ++i)
      bsSet(bs, i);
   bsClearAll(bs);
   ASSERT_EQ_INT(bs->capacity, 0);
   for (size_t i = 0; i < SIZE; ++i)
      ASSERT_EQ_INT(bsTest(bs, i), 0);
   bsFree(bs);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Fixed_ReadWrite()
{
   static size_t SIZE = 100;
   static size_t WINDOW = 3;

   BitString_t *bs = bsAlloc(SIZE);
   for (size_t wndEnd = 0; wndEnd <= SIZE + WINDOW; ++wndEnd)
   {
      // check all bits
      for (size_t i = 0; i < SIZE; ++i) {
         int value = bsTest(bs, i);
         int expectedValue = (i + WINDOW >= wndEnd && i < wndEnd) ? 1 : 0;
         //printf("wnd=%d..%d  bit[%u] act=%d exp=%d\n",
         //      (int)wndEnd - (int)WINDOW, (int)wndEnd - 1,
         //      (unsigned)i, value, expectedValue);
         ASSERT_EQ_INT(value, expectedValue);
      }

      // shift window
      if (wndEnd < SIZE)
         bsSet(bs, wndEnd);
      if (wndEnd >= WINDOW && wndEnd < SIZE + WINDOW)
         bsClear(bs, wndEnd - WINDOW);
   }
   bsFree(bs);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Variable_ReadWrite()
{
   static size_t SIZE = 100;
   static size_t WINDOW = 3;

   BitString_t *bs = bsAllocEmpty();
   for (size_t wndEnd = 0; wndEnd <= SIZE + WINDOW; ++wndEnd)
   {
      // Check all bits
      for (size_t i = 0; i < SIZE; ++i) {
         int value = bsTest(bs, i);
         int expectedValue = (i + WINDOW >= wndEnd && i < wndEnd) ? 1 : 0;
         //printf("wnd=%d..%d  bit[%u] act=%d exp=%d\n",
         //      (int)wndEnd - (int)WINDOW, (int)wndEnd - 1, (unsigned)i, value, expectedValue);
         ASSERT_EQ_INT(value, expectedValue);
      }

      // shift window
      if (wndEnd < SIZE)
         bsSet(bs, wndEnd);
      if (wndEnd >= WINDOW && wndEnd < SIZE + WINDOW)
         bsClear(bs, wndEnd - WINDOW);
   }
   bsFree(bs);
   return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Variable_Trim()
{
   BitString_t* bs = bsAllocEmpty();
   size_t size = 300;
   bsSet(bs, size - 1);
   for (; size > 0; --size) {
      ASSERT(bs->capacity < size + sizeof(long) * 8);
      if (size > 1)
         bsSet(bs, size - 2);
      bsClear(bs, size - 1);
      bsTrim(bs);
   }
   bsFree(bs);
   return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Fixed_Trim()
{
   BitString_t* bs = bsAlloc(100);
   bsTrim(bs);
   ASSERT_EQ_INT(bs->size, 100);
   bsFree(bs);
   return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Fixed_BitwiseOps()
{
   const size_t SIZE = 300;
   BitString_t* aAnd = bsAlloc(SIZE);
   BitString_t* aOr = bsAlloc(SIZE);
   BitString_t* aMinus = bsAlloc(SIZE);
   for (size_t i = 0; i < SIZE; i += 3) {
      bsSet(aAnd, i);
      bsSet(aOr, i);
      bsSet(aMinus, i);
   }
   BitString_t* b = bsAlloc(SIZE);
   for (size_t i = 0; i < SIZE; i += 5)
      bsSet(b, i);

   bsAnd(aAnd, b);
   bsOr(aOr, b);
   bsMinus(aMinus, b);

   for (size_t i = 0; i < SIZE; ++i) {
      ASSERT_EQ_INT(bsTest(aAnd, i), (i % 3 == 0) && (i % 5 == 0) ? 1 : 0);
      ASSERT_EQ_INT(bsTest(aOr, i), (i % 3 == 0) || (i % 5 == 0) ? 1 : 0);
      ASSERT_EQ_INT(bsTest(aMinus, i), (i % 3 == 0) && (i % 5 != 0) ? 1 : 0);
   }

   bsFree(b);
   bsFree(aMinus);
   bsFree(aOr);
   bsFree(aAnd);
   return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Fixed_BitwiseOps_AbortsOnSizeMismatch()
{
   const size_t SIZE = 300;
   BitString_t* a = bsAlloc(SIZE);
   BitString_t* b = bsAlloc(SIZE+1);
   ASSERT_ABORT(bsAnd(a, b));
   ASSERT_ABORT(bsOr(a, b));
   ASSERT_ABORT(bsMinus(a, b));
   return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_BitwiseOps_AbortsOnTypeMismatch()
{
   BitString_t* fixed = bsAlloc(10);
   BitString_t* variable = bsAllocEmpty();
   ASSERT_ABORT(bsAnd(fixed, variable));
   ASSERT_ABORT(bsAnd(variable, fixed));
   ASSERT_ABORT(bsOr(fixed, variable));
   ASSERT_ABORT(bsOr(variable, fixed));
   ASSERT_ABORT(bsMinus(fixed, variable));
   ASSERT_ABORT(bsMinus(variable, fixed));
   return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Variable_BitwiseAnd()
{
   const size_t SIZE = 100;
   const size_t LARGE_SIZE = SIZE + SIZE/3;
   const size_t SMALL_SIZE = SIZE - SIZE/3;
   BitString_t* largeDest = bsAllocEmpty();
   for (size_t i = 0; i < LARGE_SIZE; i += 3)
      bsSet(largeDest, i);
   BitString_t* smallDest = bsAllocEmpty();
   for (size_t i = 0; i < SMALL_SIZE; i += 3)
      bsSet(smallDest, i);

   BitString_t* src = bsAllocEmpty();
   for (size_t i = 0; i < SIZE; i += 5)
      bsSet(src, i);

   bsAnd(largeDest, src);
   bsAnd(smallDest, src);

   // Destination is larger: missing bits in src are teated as 0.
   for (size_t i = 0; i < SIZE; ++i) {
      ASSERT_EQ_INT(bsTest(largeDest, i), (i % 3 == 0) && (i % 5 == 0) ? 1 : 0);
   }
   for (size_t i = SIZE; i < LARGE_SIZE; ++i) {
      ASSERT_EQ_INT(bsTest(largeDest, i), 0);
   }

   // Destination is smaller: extra bits in src are ignored.
   for (size_t i = 0; i < SMALL_SIZE; ++i) {
      ASSERT_EQ_INT(bsTest(smallDest, i), (i % 3 == 0) && (i % 5 == 0) ? 1 : 0);
   }
   for (size_t i = SMALL_SIZE; i < SIZE; ++i) {
      ASSERT_EQ_INT(bsTest(smallDest, i), 0);
   }

   bsFree(src);
   bsFree(largeDest);
   bsFree(smallDest);
   return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

static int testCompare(BitString_t* a, BitString_t* b, size_t size)
{
   for (size_t pos = 0; pos + 1 < size; ++pos) {
      ASSERT_EQ_INT(bsCompare(a,b), 0);    // ...00    == ...00
      bsSet(a, pos);
      ASSERT_EQ_INT(bsCompare(a,b), 1);    // ...10    >  ...00
      bsSet(b, pos);
      ASSERT_EQ_INT(bsCompare(a,b), 0);    // ...10    == ...10
      bsClear(a, pos);
      ASSERT_EQ_INT(bsCompare(a,b), -1);   // ...00    <  ...10
      bsSet(a, pos + 1);
      ASSERT_EQ_INT(bsCompare(a,b), -1);   // ...01    <  ...10
      bsClear(b, pos);
      ASSERT_EQ_INT(bsCompare(a,b), 1);    // ...01    >  ...00
      bsClear(a, pos + 1);                 // ...00    == ...00
   }
   return 0;
}

TstResult BitString_Fixed_Compare()
{
   static const int SIZES[] = { 10, 100, 1000, -1 };
   int result = 0;
   for (const int* size = SIZES; result == 0 && *size > 0; ++size) {
      BitString_t* a = bsAlloc(*size);
      BitString_t* b = bsAlloc(*size);
      result |= testCompare(a, b, *size);
      if (result != 0) {
         bsPrint("a", a);
         bsPrint("b", b);
      }
      bsFree(b);
      bsFree(a);
   }
   return result;
}

TstResult BitString_Variable_Compare()
{
   static const int SIZES[] = { 10, 100, 1000, -1 };
   int result = 0;
   for (const int* size = SIZES; result == 0 && *size > 0; ++size) {
      BitString_t* a = bsAllocEmpty();
      BitString_t* b = bsAllocEmpty();
      result |= testCompare(a, b, *size);
      if (result != 0) {
         bsPrint("a", a);
         bsPrint("b", b);
      }
      bsFree(b);
      bsFree(a);
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void randomize(BitString_t* bs, size_t size)
{
   for (size_t i = 0; i < size; ++i) {
      if (mtxRandomInt(2) != 0)
         bsSet(bs, i);
      else
         bsClear(bs, i);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int testCopy(BitString_t *a, BitString_t *b, size_t size)
{
   randomize(a, size);
   randomize(b, size);
   bsCopy(b,a);
   ASSERT_EQ_INT(bsCompare(a,b), 0);
   BitString_t *c = bsDup(a);
   ASSERT_EQ_INT(bsCompare(a,c), 0);
   bsFree(c);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Fixed_Copy()
{
   const int size = 500;
   BitString_t *a, *b;
   a = bsAlloc(size);
   b = bsAlloc(size);
   int result = testCopy(a,b,size);
   bsFree(b);
   bsFree(a);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Variable_Copy()
{
   BitString_t *a, *b;
   a = bsAllocEmpty();
   b = bsAllocEmpty();
   int result = testCopy(a,b,500);
   bsFree(b);
   bsFree(a);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int testFileIo(int variable)
{
   const char FILE_NAME[] = "check.1";
   MtxFile_t* file = mfOpen(FILE_NAME,"w+b");
   const int SIZES[] = {0, 1, 2, 10, 100, 1000, 10000, -1};
   
   RngReset();
   for (const int* size = SIZES; *size > 0; ++size) {
      BitString_t* bs = variable ? bsAllocEmpty() : bsAlloc(*size);
      randomize(bs, *size);
      bsWrite(bs, file);
      bsFree(bs);
   }

   sysFseek(file->file, 0);
   RngReset();
   for (const int* size = SIZES; *size > 0; ++size) {
      BitString_t* expectedBs = variable ? bsAllocEmpty() : bsAlloc(*size);
      randomize(expectedBs, *size);
      BitString_t *bs = bsRead(file);
      ASSERT_EQ_INT(bsCompare(bs, expectedBs), 0);
      bsFree(bs);
      bsFree(expectedBs);
   }
   mfClose(file);
   sysRemoveFile(FILE_NAME);
   return 0;
}

TstResult BitString_Fixed_FileIo()
{
   return testFileIo(0);
}

TstResult BitString_Variable_FileIo()
{
   return testFileIo(1);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

static int testIntersectionCount(BitString_t *a, BitString_t *b, size_t size)
{
   size_t expectedCount = 0;
   for (size_t i = 0; i < size; ++i) {
      if (bsTest(a,i) && bsTest(b,i)) {
         ++expectedCount;
      }
   }
   ASSERT_EQ_INT(bsIntersectionCount(a,b), expectedCount);
   return 0;
}

TstResult BitString_Fixed_IntersectionCount()
{
   int result = 0;
   for (int i = 0; result == 0 && i < 10; ++i) {
      const int size = mtxRandomInt(200) + 100;
      BitString_t *a = bsAlloc(size);
      randomize(a, size);
      BitString_t *b = bsAlloc(size);
      randomize(b, size);
      result |= testIntersectionCount(a,b, size);
      bsFree(a);
      bsFree(b);
   }
   return result;
}

TstResult BitString_Variable_IntersectionCount()
{
   int result = 0;
   for (int i = 0; result == 0 && i < 10; ++i) {
      int sizeA = mtxRandomInt(200);
      int sizeB = mtxRandomInt(200);
      BitString_t *a = bsAllocEmpty();
      randomize(a, sizeA);
      BitString_t *b = bsAllocEmpty();
      randomize(b, sizeB);
      result |= testIntersectionCount(a, b, sizeA > sizeB ? sizeA : sizeB);
      bsFree(a);
      bsFree(b);
   }
   return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult BitString_Fixed_IsSubset()
{
   const size_t SIZE = 200;
   BitString_t *a = bsAlloc(SIZE);
   BitString_t *b = bsAlloc(SIZE);
   
   ASSERT(bsIsSub(a,b));
   bsSet(b, 0);
   ASSERT(bsIsSub(a,b));
   bsSet(a, 1);
   ASSERT(!bsIsSub(a,b));
   bsSet(b, 1);
   ASSERT(bsIsSub(a,b));
   bsFree(a);
   bsFree(b);
   return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

static int testIterate(BitString_t* bs, const size_t size)
{
   for (size_t start = 0; start < size / 8; ++start) {
      RngReset();
      bsClearAll(bs);
      for (int i = start; i < size; i += (RngNext() % 5 + 1)) {
         bsSet(bs, i);
      }
      size_t pos = 42;
      size_t i = start;
      ASSERT(bsFirst(bs, &pos) == 1);
      ASSERT_EQ_INT(pos, i);
      RngReset();
      while ((i += RngNext() % 5 + 1) < size) {
         ASSERT(bsNext(bs, &pos) == 1);
         ASSERT_EQ_INT(pos, i);
      }
      ASSERT(bsNext(bs, &pos) == 0);
   }
   return 0;
}

TstResult BitString_Fixed_Iterate()
{
   const size_t SIZE = 400;
   BitString_t *bs = bsAlloc(SIZE);
   int result = testIterate(bs, SIZE);
   bsFree(bs);
   return result;
}

TstResult BitString_Variable_Iterate()
{
   const size_t SIZE = 400;
   BitString_t *bs = bsAllocEmpty();
   int result = testIterate(bs, SIZE);
   bsFree(bs);
   return result;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
