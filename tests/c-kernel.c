////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for kernel functions
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"

#include <stdlib.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestField1()
{
   FEL a,b,c;
   int ai, bi, ci;

   /* Test one and zero elements
      -------------------------- */
   for (ai = 0; ai < (int)FfOrder; ++ai) {
      a = FTab[ai];
      if (FfAdd(a,FF_ZERO) != a) {
         TST_FAIL2("%d+0=%d\n",(int)a,(int)FfAdd(a,FF_ZERO));
      }
      if (FfMul(a,FF_ONE) != a) {
         TST_FAIL2("%d*1=%d\n",(int)a,(int)FfMul(a,FF_ONE));
      }
   }

   /* Test negative and inverse
      ------------------------- */
   for (ai = 0; ai < (int)FfOrder; ++ai) {
      a = FTab[ai];
      b = FfNeg(a);
      if (!ISFEL(b) || (FfAdd(a,b) != FF_ZERO)) {
         TST_FAIL2("Illegal negative: -%d = %d",a,b);
      }
      if (a != FF_ZERO) {
         b = FfInv(a);
         if (!ISFEL(b) || (FfMul(a,b) != FF_ONE)) {
            TST_FAIL("Illegal inverse");
         }
      }
   }

   /* Test arithmetic
      --------------- */
   for (ai = 0; ai < (int)FfOrder; ++ai) {
      a = FTab[ai];
      for (bi = ai; bi < (int)FfOrder; ++bi) {
         b = FTab[bi];
         c = FfAdd(a,b);
         if (!ISFEL(c)) { TST_FAIL("FfAdd() error"); }
         if (c != FfAdd(b,a)) { TST_FAIL("'+' not commutative"); }
         c = FfMul(a,b);
         if (!ISFEL(c)) { TST_FAIL("FfMul() error"); }
         if (c != FfMul(b,a)) { TST_FAIL("'*' not commutative"); }

         for (ci = 0; ci < (int)FfOrder; ++ci) {
            c = FTab[ci];
            if (FfAdd(a,FfAdd(b,c)) != FfAdd(FfAdd(a,b),c)) {
               TST_FAIL("'+' not associative");
            }
            if (FfMul(a,FfMul(b,c)) != FfMul(FfMul(a,b),c)) {
               TST_FAIL("'*' not associative");
            }
            if (FfMul(a,FfAdd(b,c)) != FfAdd(FfMul(a,b),FfMul(a,c))) {
               TST_FAIL("a*(b+c) != a*b+a*c");
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F FiniteFieldArithmetic()
{
   while (NextField() > 0) {
      TestField1();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestGen1()
{
   FEL a, b;
   int i;

   a = FfGen;
   b = a;
   for (i = 1; i < (int)FfOrder - 1; ++i) {
      if (b == FF_ONE) {
         TST_FAIL1("Generator has order %d",i);
      }
      b = FfMul(a,b);
   }
   if (b != FF_ONE) {
      TST_FAIL2("g^%d = %d, should be %d",b,FF_ONE);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F FfGenerator()
{
   while (NextField() > 0) {
      TestGen1();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestInsertExtract2(PTR x, int pos)
{
   int i1, i2, i3;
   int max1 = (FfOrder > 32) ? 32 : FfOrder;
   int max3 = (FfOrder > 16) ? 16 : FfOrder;

   for (i1 = 0; i1 < max1; ++i1) {
      FfInsert(x,pos + 1 - 1,FTab[i1]);
      for (i2 = 0; i2 < FfOrder; ++i2) {
         FfInsert(x,pos + 2 - 1,FTab[i2]);
         for (i3 = 0; i3 < max3; ++i3) {
            FfInsert(x,pos + 3 - 1,FTab[i3]);

            if ((FfExtract(x,pos + 1 - 1) != FTab[i1]) ||
                (FfExtract(x,pos + 2 - 1) != FTab[i2]) ||
                (FfExtract(x,pos + 3 - 1) != FTab[i3])
                ) {
               TST_FAIL6("Read %d %d %d, expected %d %d %d",
                     FfExtract(x,pos + 1 - 1),FfExtract(x,pos + 2 - 1),
                     FfExtract(x,pos + 3 - 1),FTab[i1],FTab[i2],FTab[i3]);
            }

         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestInsertExtract1()
{
   int pos;
   PTR x;

   FfSetNoc(20);
   x = FfAlloc(1);
   for (pos = 0; pos < 14; ++pos) {
      TestInsertExtract2(x,pos);
   }
   free(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F InsertExtract()
{
   while (NextField() > 0) {
      TestInsertExtract1();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestFindPiv2(PTR x, int noc)
{
   FEL a;
   int pos;
   int i;

   if ((pos = FfFindPivot(x,&a)) >= 0) {
      TST_FAIL2("Found pivot %d at column %d in empty row",a,pos);
   }
   for (i = 0; i < noc; ++i) {
      FfInsert(x,i,FTab[i % (FfOrder - 1) + 1]);
   }
   for (i = 0; i < noc; ++i) {
      int k;
      for (k = 1; k < FfOrder; ++k) {
         FfInsert(x,i,FTab[k]);
         if (((pos = FfFindPivot(x,&a)) != i) || (a != FTab[k])) {
            TST_FAIL4("Found Pivot %d at %d, expected %d at %d",
                  a,pos,FTab[k],i);
         }
      }
      FfInsert(x,i,FF_ZERO);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestFindPiv3(PTR x, int noc)
{
   FEL a;
   int pos;
   int i;

   FfMulRow(x,FF_ZERO);
   for (i = noc - 1; i > 0; --i) {
      FfInsert(x,i,FF_ONE);
      pos = FfFindPivot(x,&a);
      if (pos != i) {
         TST_FAIL2("Found pivot at %d, expected %d",pos,i);
      }
      FfSetNoc(i);
      pos = FfFindPivot(x,&a);
      if (pos != -1) {
         TST_FAIL2("Found pivot at %d, expected %d",pos,-1);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestFindPiv1()
{
   int noc;
   PTR x;
   for (noc = 0; noc < 35; ++noc) {
      FfSetNoc(noc);
      x = FfAlloc(1);
      TestFindPiv2(x,noc);
      TestFindPiv3(x,noc);
      SysFree(x);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F FindPivot()
{
   while (NextField() > 0) {
      TestFindPiv1();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestFelToInt1()
{
   long l, m;
   FEL a;

   ASSERT_EQ_INT(FfFromInt(0), FF_ZERO);
   ASSERT_EQ_INT(FfFromInt(1), FF_ONE);
   ASSERT_EQ_INT(FfToInt(FF_ZERO), 0);
   ASSERT_EQ_INT(FfToInt(FF_ONE), 1);
   for (l = 0; l < FfOrder; ++l) {
      FTab[l] = FfFromInt(l);
      ASSERT(ISFEL(FTab[l]));
      for (m = 0; m < l; ++m) {
         if (FTab[m] == FTab[l]) {
            TST_FAIL2("FfFromInt(%d) = zitof(%d)",m,l);
         }
      }
      ASSERT_EQ_INT(FfToInt(FTab[l]),l);
   }

   // check that char(F) = 0
   for (a = FF_ZERO, l = 0; l < FfChar; ++l, a = FfAdd(a,FF_ONE)) {
      ASSERT_EQ_INT(FTab[l], a);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F FfToInteger()
{
   while (NextField() > 0) {
      TestFelToInt1();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestSubfield1(int fld, int sub)
{
   FEL tabsub[256];
   FEL tabemb[256];
   int i;

   FfSetField(sub);
   for (i = 0; i < sub; ++i) {
      tabsub[i] = FfFromInt(i);
   }
   FfSetField(fld);
   for (i = 0; i < sub; ++i) {
      tabemb[i] = FfEmbed(tabsub[i],sub);
      ASSERT(ISFEL(tabemb[i]));
      ASSERT_EQ_INT(FfRestrict(tabemb[i],sub),tabsub[i]);
   }

   for (i = 0; i < sub; ++i) {
      int k;
      for (k = 0; k < sub; ++k) {
         int l;
         FEL x = FfAdd(tabemb[i],tabemb[k]);
         FEL y = FfMul(tabemb[i],tabemb[k]);
         for (l = 0; l < sub && tabemb[l] != x; ++l) {
         }
         if (l >= sub) {
            TST_FAIL2("Embedding of GF(%d) into GF(%d) not closed (+)",fld,sub);
         }
         for (l = 0; l < sub && tabemb[l] != y; ++l) {
         }
         if (l >= sub) {
            TST_FAIL2("Embedding of GF(%d) into GF(%d) not closed (*)",fld,sub);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F Subfields()
{
   TestSubfield1(256,2);
   TestSubfield1(256,4);
   TestSubfield1(256,16);
   TestSubfield1(128,2);
   TestSubfield1(64,2);
   TestSubfield1(64,4);
   TestSubfield1(64,8);
   TestSubfield1(32,2);
   TestSubfield1(16,2);
   TestSubfield1(16,4);

   TestSubfield1(243,3);
   TestSubfield1(81,3);
   TestSubfield1(81,9);
   TestSubfield1(9,3);

   TestSubfield1(125,5);
   TestSubfield1(25,5);

   TestSubfield1(49,7);

   TestSubfield1(121,11);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestAddRow2(PTR x, PTR y, int noc, int d1, int d2)
{
   int i, k;

   for (i = 0; i < noc; ++i) {
      FfInsert(x,i,FTab[(i + d1) % FfOrder]);
      FfInsert(y,i,FTab[(i + d2) % FfOrder]);
   }
   FfAddRow(x,y);
   for (i = 0; i < noc; ++i) {
      FEL f = FfExtract(x,i);
      if (f != FfAdd(FTab[(i + d1) % FfOrder],FTab[(i + d2) % FfOrder])) {
         TST_FAIL3("FfAddRow failed at %d+%d in column %d",
               FTab[(i + d1) % FfOrder],FTab[(i + d2) % FfOrder],i);
      }
   }

   for (k = 0; k < FfOrder; ++k) {
      for (i = 0; i < noc; ++i) {
         FfInsert(x,i,FTab[(i + d1) % FfOrder]);
      }
      FfAddMulRow(x,y,FTab[k]);
      for (i = 0; i < noc; ++i) {
         FEL f = FfExtract(x,i);
         FEL g = FfAdd(FTab[(i + d1) % FfOrder],
                       FfMul(FTab[(i + d2) % FfOrder],FTab[k]));
         if (f != g) {
            TST_FAIL4("FfAddMulRow failed at %d+%d*%d in column %d",
                  FTab[(i + d1) % FfOrder],FTab[k],FTab[(i + d2) % FfOrder],i);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestAddRow1a(PTR x, PTR y, int noc)
{
   int i;
   for (i = 0; i < (noc > 8 ? 8 : noc); ++i) {
      int ai;
      for (ai = 0; ai < FfOrder; ++ai) {
         int bi;
         FfInsert(x,i,FTab[ai]);
         for (bi = 0; bi < FfOrder; ++bi) {
            FEL ist, soll;
            FfInsert(y,i,FTab[bi]);
            FfAddRow(y,x);
            ist = FfExtract(y,i);
            soll = FfAdd(FTab[ai],FTab[bi]);
            if (ist != soll) {
               TST_FAIL5("FfAddrow failed at col %d: %d+%d=%d, expected %d",
                     i,FTab[ai],FTab[bi],ist,soll);
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestAddRow1()
{
   const int noc = 16;
   int max1 = FfOrder < 32 ? FfOrder : 32;
   PTR x, y;
   int i, k;

   FfSetNoc(noc);
   x = FfAlloc(1);
   y = FfAlloc(1);
   TestAddRow1a(x,y,noc);
   for (i = 0; i < max1; ++i) {
      for (k = 0; k < FfOrder; ++k) {
         TestAddRow2(x,y,noc,i,k);
      }
   }
   free(x);
   free(y);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F RowOperations()
{
   while (NextField() > 0) {
      TestAddRow1();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMulRow1(PTR x, int noc)
{
   int i;
   FEL a;

   for (a = FF_ONE, i = 1; i < FfOrder; a = FfMul(a,FTab[i++])) {
   }
   a = FfInv(a);

   for (i = 0; i < 500; ++i) {
      int val, l;
      for (val = i, l = 0; l < noc; ++l) {
         FfInsert(x,l,FTab[val % FfOrder]);
         if ((val /= FfOrder) == 0) {
            val = i;
         }
      }
      for (l = 1; l < FfOrder; ++l) {
         FfMulRow(x,FTab[l]);
      }
      FfMulRow(x,a);
      for (val = i, l = 0; l < noc; ++l) {
         if (FfExtract(x,l) != FTab[val % FfOrder]) {
            TST_FAIL("FfMulRow() failed");
         }
         if ((val /= FfOrder) == 0) {
            val = i;
         }
      }
   }

   FfMulRow(x,FF_ZERO);
   for (i = 0; i < noc; ++i) {
      if (FfExtract(x,i) != FF_ZERO) {
         TST_FAIL("FfMulRow(x,0) failed");
      }
   }

}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F MulRow()
{
   while (NextField() > 0) {
      int noc;
      for (noc = 0; noc < 35; ++noc) {
         PTR x;
         FfSetNoc(noc);
         x = FfAlloc(1);
         TestMulRow1(x,noc);
         free(x);
      }
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestPtr2(PTR x)
{
   PTR xx;
   xx = x;
   FfStepPtr(&xx);
   if ((char *)xx != (char *)x + FfCurrentRowSize) {
      TST_FAIL("Pointer increment error");
   }
}


static void TestPtr1()
{
   int noc;
   PTR x;

   for (noc = 10; noc < 30; ++noc) {
      FfSetNoc(noc);
      x = FfAlloc(noc);
      TestPtr2(x);
      free(x);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F RowPointerOperations()
{
   while (NextField() > 0) {
      TestPtr1();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestRowSize1()
{
   int noc;
   for (noc = 0; noc < 100; ++noc) {
      int rs;
      int diff;
      rs = FfRowSize(noc);
      if ((rs < 0) || (rs > noc + sizeof(long))) {
         TST_FAIL2("FfRowSize(%d) = %d out of range",noc,rs);
      }
      diff = rs - FfTrueRowSize(noc);
      if ((diff < 0) || (diff >= sizeof(long))) {
         TST_FAIL("FfRowSize() and FfTrueRowSize() differ too much");
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F RowSize()
{
   while (NextField() > 0) {
      TestRowSize1();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestCmpRows2(PTR m1, PTR m2, int noc)
{
   int i;

   FfMulRow(m1,FF_ZERO);
   FfMulRow(m2,FF_ZERO);

   for (i = 1; i < FfOrder; ++i) {
      int k;
      for (k = 0; k < noc; ++k) {
         if (FfCmpRows(m2,m1) != 0) {
            TST_FAIL("Rows are different");
         }
         FfInsert(m1,k,FTab[i]);
         if (FfCmpRows(m2,m1) == 0) {
            TST_FAIL("Rows are still equal");
         }
         FfInsert(m2,k,FTab[i]);
         if (FfCmpRows(m2,m1) != 0) {
            TST_FAIL("Rows are still different");
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestCmpRows1()
{
   int noc;

   for (noc = 10; noc < 30; ++noc) {
      PTR m1, m2;
      FfSetNoc(noc);
      m1 = FfAlloc(1);
      m2 = FfAlloc(1);
      TestCmpRows2(m1,m2,noc);
      free(m1);
      free(m2);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F CompareRows(unsigned flags)
{
   while (NextField() > 0) {
      TestCmpRows1();
   }
}
