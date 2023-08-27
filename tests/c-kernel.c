////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for kernel functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>


////////////////////////////////////////////////////////////////////////////////////////////////////

static int checkIntMapping(char* isSet, int bufSize)
{
   ASSERT_EQ_INT(ffFromInt(0), FF_ZERO);
   ASSERT_EQ_INT(ffFromInt(1), FF_ONE);
   ASSERT_EQ_INT(ffToInt(FF_ZERO), 0);
   ASSERT_EQ_INT(ffToInt(FF_ONE), 1);

   // Check that ffFromInt() and ffToInt() are inverse
   memset(isSet, 0, bufSize);
   for (int i = 0; i < ffOrder; ++i) {
      FEL f = ffFromInt(i);
      ASSERT(f < bufSize);
      ASSERT(ISFEL(f));
      ASSERT_EQ_INT(isSet[f], 0);
      isSet[f] = 1;
      ASSERT_EQ_INT(ffToInt(f), i);
   }
   return 0;
}

TstResult Kernel_Field_IntFelMapping(int q)
{
   char* isSet = NALLOC(char, 0x10000);
   int result = checkIntMapping(isSet, 0x10000);
   sysFree(isSet);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_Field_Characteristic(int q)
{
   FEL a = FF_ZERO;
   for (int i = 1; i <= ffChar; ++i) {
      a = ffAdd(a, FF_ONE);
      if (i < ffChar)
         ASSERT(a != FF_ZERO);
      else 
         ASSERT(a == FF_ZERO);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_Field_Inversion(int q)
{
   // one and zero
   for (int ai = 0; ai < (int)ffOrder; ++ai) {
      const FEL a = FTab[ai];
      ASSERT_EQ_INT(ffAdd(a,FF_ZERO),a);
      ASSERT_EQ_INT(ffMul(a,FF_ONE),a);
   }

   // additive and multiplicative inverse
   for (int ai = 0; ai < (int)ffOrder; ++ai) {
      const FEL a = FTab[ai];
      const FEL minusA = ffNeg(a);
      ASSERT(ISFEL(minusA));
      ASSERT_EQ_INT(ffAdd(a,minusA),FF_ZERO);

      if (a != FF_ZERO) {
	 const FEL inverseOfA = ffInv(a);
         ASSERT(ISFEL(inverseOfA));
	 ASSERT_EQ_INT(ffMul(a,inverseOfA), FF_ONE);

	 const FEL b = ffAdd(a, FF_ONE);
	 const FEL x = ffDiv(ffAdd(a, FF_ONE), a);
	 const FEL y = ffMul(a, x);
	 ASSERT_EQ_INT(y, b);
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_Field_Distributivity(int q)
{
   //puts(__func__);
   const int step = ffOrder / 257 + 1;

   for (int ai = 0; ai < (int)ffOrder; ai += step) {
      const FEL a = FTab[ai];
      for (int bi = ai; bi < (int)ffOrder; bi += step) {
         const FEL b = FTab[bi];
         for (int ci = ci; ci < (int)ffOrder; ci += step) {
            const FEL c = FTab[ci];

            const FEL ab = ffMul(a,b);
            const FEL ac = ffMul(a,c);
            ASSERT_EQ_INT(ffAdd(ab,ac), ffMul(a, ffAdd(b, c))); // a * b + a * c = a * (b + c)
         }
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_Field_CommutativityAndAssociativity(int q)
{
   const int step = ffOrder / 257 + 1;

   for (int ai = 0; ai < (int)ffOrder; ai += step) {
      const FEL a = FTab[ai];
      for (int bi = ai; bi < (int)ffOrder; bi += step) {
         const FEL b = FTab[bi];

         {
            const FEL ab = ffAdd(a,b);
            const FEL ba = ffAdd(b,a);
            ASSERT_EQ_INT(ab, ba);      // a + b = b + a
            const FEL ab_a = ffAdd(ab,a);
            const FEL a_ba = ffAdd(a,ba);
            ASSERT_EQ_INT(ab_a, a_ba);  // (a + b) + a = a + (b + a)
         }

         {
            const FEL ab = ffMul(a,b);
            const FEL ba = ffMul(b,a);
            ASSERT_EQ_INT(ab, ba);      // a * b = b * a
            const FEL ab_a = ffMul(ab,a);
            const FEL a_ba = ffMul(a,ba);
            ASSERT_EQ_INT(ab_a, a_ba);  // (a * b) * a = a * (b * a)
         }
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_Field_AddSub(int q)
{
   RngReset();
   for (int i = 0; i < ffOrder; ++i) {
      const FEL a = RandomFieldElement();
      const FEL b = RandomFieldElement();
      {
         FEL ab = ffAdd(a,b);
         FEL ab_a = ffSub(ab,a);
         ASSERT_EQ_INT(ab_a, b); // (a+b) - a = b
         FEL ab_b = ffSub(ab,b);
         ASSERT_EQ_INT(ab_b, a); // (a+b) - b = a
      }
      {
         FEL ab = ffSub(a,b);
         FEL ba = ffSub(b,a);
         FEL ab_ba = ffAdd(ab,ba);
         ASSERT_EQ_INT(ab_ba, FF_ZERO); // (a-b) + (b-a) = 0
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_Field_MulDiv(int q)
{
   //puts(__func__);
   RngReset();
   for (int i = 0; i < ffOrder; ++i) {
      const FEL a = RandomNonzeroFieldElement();
      const FEL b = RandomNonzeroFieldElement();
      {
         FEL ab = ffMul(a,b);
         FEL ab_a = ffDiv(ab,a);
         ASSERT_EQ_INT(ab_a, b); // (a * b) / a = b
         FEL ab_b = ffDiv(ab,b);
         ASSERT_EQ_INT(ab_b, a); // (a * b) / b = a
      }
      {
         FEL ab = ffDiv(a,b);
         FEL ba = ffDiv(b,a);
         FEL ab_ba = ffMul(ab,ba);
         ASSERT_EQ_INT(ab_ba, FF_ONE); // (a/b) * (b/a) = 0
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_Field_Generator(int q)
{
   FEL b = ffGen;
   for (int i = 1; i < ffOrder - 1; ++i) {
      ASSERT(b != FF_ONE);
      b = ffMul(b,ffGen);
   }
   ASSERT_EQ_INT(b, FF_ONE);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int InsertExtract2(PTR x, int pos)
{
   const int max1 = (ffOrder > 32) ? 32 : ffOrder;
   const int max3 = (ffOrder > 16) ? 16 : ffOrder;

   for (int i1 = 0; i1 < max1; ++i1) {
      ffInsert(x,pos + 1 - 1,FTab[i1]);
      for (int i2 = 0; i2 < ffOrder; ++i2) {
         ffInsert(x,pos + 2 - 1,FTab[i2]);
         for (int i3 = 0; i3 < max3; ++i3) {
            ffInsert(x,pos + 3 - 1,FTab[i3]);
	    ASSERT_EQ_INT(ffExtract(x,pos + 0), FTab[i1]);
	    ASSERT_EQ_INT(ffExtract(x,pos + 1), FTab[i2]);
	    ASSERT_EQ_INT(ffExtract(x,pos + 2), FTab[i3]);
         }
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_RowOps_InsertExtract(int q)
{
   int result = 0;
   const PTR x = ffAlloc(1, 20);
   for (int pos = 0; result == 0 && pos < 14; ++pos) {
      result |= InsertExtract2(x,pos);
   }
   free(x);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestFindPiv2(PTR row, int noc)
{
   // fill with nonzero elemets
   for (int i = 0; i < noc; ++i) {
      ffInsert(row,i,FTab[i % (ffOrder - 1) + 1]);
   }

   // test each column
   for (int i = 0; i < noc; ++i) {
      for (int k = 1; k < ffOrder; ++k) {
         ffInsert(row,i,FTab[k]);
	 FEL pivotElement;
	 const int pivotColumn = ffFindPivot(row,&pivotElement, noc);
         ASSERT_EQ_INT(pivotColumn, i);
	 ASSERT_EQ_INT(pivotElement, FTab[k]);
      }
      ffInsert(row,i,FF_ZERO);
   }

   // empty row
   FEL dummy;
   ASSERT_EQ_INT(ffFindPivot(row,&dummy, noc), -1);
   return 0;
}

static int TestFindPiv3(PTR row, int noc)
{
   int i;

   ffMulRow(row,FF_ZERO, noc);
   for (i = noc - 1; i > 0; --i) {
      ffInsert(row,i,FF_ONE);
      FEL pivotValue;
      const int pivotColumn = ffFindPivot(row,&pivotValue, noc);
      ASSERT_EQ_INT(pivotColumn , i);

      // reduce row size below pivot column and try again
      ASSERT_EQ_INT(ffFindPivot(row,&pivotValue, i), -1);
   }
   return 0;
}

TstResult Kernel_FindPivot(int q)
{
   int result = 0;
   for (int noc = 0; result == 0 && noc < 35; ++noc) {
      const PTR x = ffAlloc(1, noc);
      result |= TestFindPiv2(x,noc);
      result |= TestFindPiv3(x,noc);
      sysFree(x);
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestSubfield1(int fld, int sub)
{
   FEL tabsub[256];
   FEL tabemb[256];

   // Fill table with all subfield elememts
   ffSetField(sub);
   for (int i = 0; i < sub; ++i) {
      tabsub[i] = ffFromInt(i);
   }

   // Switch to main field and calculate embedded subfield elements.
   // Verify that embedding and restriction are compatible.
   ffSetField(fld);
   for (int i = 0; i < sub; ++i) {
      tabemb[i] = ffEmbed(tabsub[i],sub);
      ASSERT(ISFEL(tabemb[i]));
      ASSERT_EQ_INT(ffRestrict(tabemb[i],sub),tabsub[i]);
   }

   // Verify that the embedded subfield is closed under both '+' and '*'.
   for (int i = 0; i < sub; ++i) {
      for (int k = 0; k < sub; ++k) {
         const FEL sum = ffAdd(tabemb[i],tabemb[k]);
         for (int l = 0; tabemb[l] != sum; ++l) {
            if (l >= sub)
               TST_FAIL("Embedding of F%d into F%d not closed (+)", sub, fld);
         }
         const FEL prod = ffMul(tabemb[i],tabemb[k]);
         for (int l = 0; tabemb[l] != prod; ++l) {
            if (l >= sub)
               TST_FAIL("Embedding of F%d into F%d not closed (*)", sub, fld);
         }
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_Field_Subfields()
{
   int result = 0;
   result |= TestSubfield1(256,2);
   result |= TestSubfield1(256,4);
   result |= TestSubfield1(256,16);
   result |= TestSubfield1(128,2);
   result |= TestSubfield1(64,2);
   result |= TestSubfield1(64,4);
   result |= TestSubfield1(64,8);
   result |= TestSubfield1(16,2);
   result |= TestSubfield1(16,4);

   result |= TestSubfield1(243,3);
   result |= TestSubfield1(81,3);
   result |= TestSubfield1(81,9);
   result |= TestSubfield1(9,3);

   result |= TestSubfield1(125,5);
   result |= TestSubfield1(25,5);

   result |= TestSubfield1(49,7);

   result |= TestSubfield1(121,11);

   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestAddRow2(PTR x, PTR y, int noc, int d1, int d2)
{
   int i, k;

   for (i = 0; i < noc; ++i) {
      ffInsert(x,i,FTab[(i + d1) % ffOrder]);
      ffInsert(y,i,FTab[(i + d2) % ffOrder]);
   }
   ffAddRow(x,y, noc);
   for (i = 0; i < noc; ++i) {
      FEL f = ffExtract(x,i);
      if (f != ffAdd(FTab[(i + d1) % ffOrder],FTab[(i + d2) % ffOrder])) {
         TST_FAIL("FfAddRow failed at %d+%d in column %d",
               FTab[(i + d1) % ffOrder],FTab[(i + d2) % ffOrder],i);
      }
   }

   const int step = (ffOrder <= 256) ? 1 : ffOrder / 100;
   for (k = 0; k < ffOrder; k += step) {
      for (i = 0; i < noc; ++i) {
         ffInsert(x,i,FTab[(i + d1) % ffOrder]);
      }
      ffAddMulRow(x,y,FTab[k],noc);
      for (i = 0; i < noc; ++i) {
         FEL f = ffExtract(x,i);
         FEL g = ffAdd(FTab[(i + d1) % ffOrder],
                       ffMul(FTab[(i + d2) % ffOrder],FTab[k]));
         if (f != g) {
            TST_FAIL("FfAddMulRow failed at %d+%d*%d in column %d",
                  FTab[(i + d1) % ffOrder],FTab[k],FTab[(i + d2) % ffOrder],i);
         }
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestAddRowPartial(PTR x, PTR y, int noc, int d1, int d2)
{
   int i, k;

   const int step = (ffOrder <= 256) ? 1 : ffOrder / 100;
   for (i = 0; i < noc; ++i) {
      ffInsert(x,i,FTab[(i + d1) % ffOrder]);
      ffInsert(y,i,FTab[(i + d2) % ffOrder]);
   }
   ffAddRow(x,y, noc);
   for (i = 0; i < noc; ++i) {
      FEL f = ffExtract(x,i);
      if (f != ffAdd(FTab[(i + d1) % ffOrder],FTab[(i + d2) % ffOrder])) {
         TST_FAIL("FfAddRow failed at %d+%d in column %d",
               FTab[(i + d1) % ffOrder],FTab[(i + d2) % ffOrder],i);
      }
   }

   for (k = 0; k < ffOrder; k += step) {
      for (i = 0; i < noc; ++i) {
         ffInsert(x,i,FTab[(i + d1) % ffOrder]);
      }
      ffAddMulRow(x,y,FTab[k], noc);
      for (i = 0; i < noc; ++i) {
         FEL f = ffExtract(x,i);
         FEL g = ffAdd(FTab[(i + d1) % ffOrder],
                       ffMul(FTab[(i + d2) % ffOrder],FTab[k]));
         if (f != g) {
            TST_FAIL("FfAddMulRow failed at %d+%d*%d in column %d",
                  FTab[(i + d1) % ffOrder],FTab[k],FTab[(i + d2) % ffOrder],i);
         }
      }
   }
   return 0;
}

static int TestAddRow1a(PTR x, PTR y, int noc)
{
   const int step = (ffOrder <= 256) ? 1 : ffOrder / 100;
   for (int i = 0; i < (noc > 8 ? 8 : noc); ++i) {
      for (int ai = 0; ai < ffOrder; ai += step) {
         ffInsert(x,i,FTab[ai]);
         for (int bi = 0; bi < ffOrder; bi += step) {
            FEL ist, soll;
            ffInsert(y,i,FTab[bi]);
            ffAddRow(y,x, noc);
            ist = ffExtract(y,i);
            soll = ffAdd(FTab[ai],FTab[bi]);
            if (ist != soll) {
               TST_FAIL("FfAddrow failed at col %d: %d+%d=%d, expected %d",
                     i,FTab[ai],FTab[bi],ist,soll);
               return 1;
            }
         }
      }
   }
   return 0;
}

TstResult Kernel_RowOps_AddRow(int q)
{
   const int noc = 16;
   const int max1 = ffOrder < 32 ? ffOrder : 32;
   int result = 0;

   PTR x = ffAlloc(1, noc);
   PTR y = ffAlloc(1, noc);
   TestAddRow1a(x,y,noc);
   for (int i = 0; result == 0 && i < max1; ++i) {
      const int step = (ffOrder <= 256) ? 1 : ffOrder / 100;
      for (int k = 0; result == 0 && k < ffOrder; k += step) {
         result |= TestAddRow2(x,y,noc,i,k);
         result |= TestAddRowPartial(x,y,noc,i,k);
      }
   }
   free(x);
   free(y);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestMulRow1(PTR row, FEL* row2, int noc)
{
   // Fill row with random elements
   for (int i = 0; i < noc; ++i) {
      row2[i] = RandomFieldElement();
      ffInsert(row, i, row2[i]);
   }

   // Multiply row with nonzero field elements
   const int step = ffOrder / 1500 + 1;
   for (int i = 1; i < ffOrder; i += step) {
      FEL a = FTab[i];
      for (int col = 0; col < noc; ++col)
         row2[col] = ffMul(row2[col], a);
      ffMulRow(row, a, noc);

      for (int col = 0; col < noc; ++col) {
         ASSERT_EQ_INT(ffExtract(row,col), row2[col]);
      }
   }

   // Multiply with zero.
   ffMulRow(row,FF_ZERO, noc);
   for (int col = 0; col < noc; ++col) {
      ASSERT_EQ_INT(ffExtract(row,col), FF_ZERO);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_RowOps_MulRow(int q)
{
   const int noc = ffOrder + 100;
   int result = 0;
   
   PTR row = ffAlloc(1, noc);
   FEL* row2 = NALLOC(FEL, noc);
   result = TestMulRow1(row, row2, noc);
   free(row);
   sysFree(row2);
  
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_RowOps_MulRowPadsWithZero(int q)
{
   ffSetField(2);
   PTR x = ffAlloc(1, 1);
   uint8_t* const begin = (uint8_t*) x;
   uint8_t* const pad = begin + ffRowSizeUsed(1);
   uint8_t* const end = begin + ffRowSize(1);

   memset(x,0xaa,ffRowSize(1));
   ffMulRow(x, FF_ZERO, 1);

   for (const uint8_t* p = pad; p < end; ++p)
      ASSERT_EQ_INT(*p, 0);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestStepPtr(PTR x, int noc)
{
   PTR xx = x;
   ffStepPtr(&xx, noc);
   if ((char *)xx != (char *)x + ffRowSize(noc)) {
      TST_FAIL("Pointer increment error",0);
   }
   return 0;
}

TstResult Kernel_RowOps_StepPtr(int q)
{
   int result = 0;
   for (int noc = 10; result == 0 && noc < 30; ++noc) {
      PTR x = ffAlloc(noc, noc);
      result |= TestStepPtr(x, noc);
      ffFree(x);
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_RowOps_ffRowSize_AbortsOnNegativeArgument(int q)
{
   ASSERT_ABORT(ffRowSize(-1));
   ASSERT_ABORT(ffSize(1,-1));
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_RowOps_ffSize_SupportsNegativeArgument(int q)
{
   ASSERT_EQ_INT(ffSize(-3, 20), -ffSize(3,20));
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_RowOps_RowSize(int q)
{
   for (int noc = 0; noc < 100; ++noc) {
      int rs = ffRowSize(noc);
      if ((rs < 0) || (rs > noc * sizeof(FEL) + sizeof(long))) {
         TST_FAIL("ffRowSize(%d) = %d out of range",noc,rs);
      }
      int diff = rs - ffRowSizeUsed(noc);
      if ((diff < 0) || (diff >= sizeof(long))) {
         TST_FAIL("ffRowSize() and ffRowSizeUsed() differ too much",0);
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestCmpRows2(PTR m1, PTR m2, int noc)
{
   int i;

   ffMulRow(m1,FF_ZERO, noc);
   ffMulRow(m2,FF_ZERO, noc);

   for (i = 1; i < ffOrder; ++i) {
      int k;
      for (k = 0; k < noc; ++k) {
         if (ffCmpRows(m2,m1, noc) != 0) {
            TST_FAIL("Rows are different",0);
         }
         ffInsert(m1,k,FTab[i]);
         if (ffCmpRows(m2,m1, noc) == 0) {
            TST_FAIL("Rows are still equal",0);
         }
         ffInsert(m2,k,FTab[i]);
         if (ffCmpRows(m2,m1, noc) != 0) {
            TST_FAIL("Rows are still different", 0);
         }
      }
   }
   return 0;
}

TstResult Kernel_RowOps_CmpRows(int q)
{
   int result = 0;

   for (int noc = 10; result == 0 && noc < 30; ++noc) {
      PTR m1, m2;
      m1 = ffAlloc(1, noc);
      m2 = ffAlloc(1, noc);
      result |= TestCmpRows2(m1,m2,noc);
      free(m1);
      free(m2);
   }
   return result;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
