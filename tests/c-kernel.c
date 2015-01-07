////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check functions for kernel.c.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"
#include "c-kernel.h"

#include <stdlib.h>
#include <stdio.h>

/* --------------------------------------------------------------------------
   Testfield() - Test finite field arithmetic
   -------------------------------------------------------------------------- */

static void TestField1()
{
   FEL a,b,c;
   int ai, bi, ci;

   /* Test one and zero elements
      -------------------------- */
   for (ai = 0; ai < (int)FfOrder; ++ai) {
      a = FTab[ai];
      if (FfAdd(a,FF_ZERO) != a) {
         Error("%d+0=%d\n",(int)a,(int)FfAdd(a,FF_ZERO));
      }
      if (FfMul(a,FF_ONE) != a) {
         Error("%d*1=%d\n",(int)a,(int)FfMul(a,FF_ONE));
      }
   }

   /* Test negative and inverse
      ------------------------- */
   for (ai = 0; ai < (int)FfOrder; ++ai) {
      a = FTab[ai];
      b = FfNeg(a);
      if (!ISFEL(b) || (FfAdd(a,b) != FF_ZERO)) {
         Error("Illegal negative: -%d = %d",a,b);
      }
      if (a != FF_ZERO) {
         b = FfInv(a);
         if (!ISFEL(b) || (FfMul(a,b) != FF_ONE)) {
            Error("Illegal inverse");
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
         if (!ISFEL(c)) { Error("FfAdd() error"); }
         if (c != FfAdd(b,a)) { Error("'+' not commutative"); }
         c = FfMul(a,b);
         if (!ISFEL(c)) { Error("FfMul() error"); }
         if (c != FfMul(b,a)) { Error("'*' not commutative"); }

         for (ci = 0; ci < (int)FfOrder; ++ci) {
            c = FTab[ci];
            if (FfAdd(a,FfAdd(b,c)) != FfAdd(FfAdd(a,b),c)) {
               Error("'+' not associative");
            }
            if (FfMul(a,FfMul(b,c)) != FfMul(FfMul(a,b),c)) {
               Error("'*' not associative");
            }
            if (FfMul(a,FfAdd(b,c)) != FfAdd(FfMul(a,b),FfMul(a,c))) {
               Error("a*(b+c) != a*b+a*c");
            }
         }
      }
   }
}


void TestField(unsigned flags)
{
   while (NextField() > 0) {
      TestField1();
   }
   flags = 0;
}


/* --------------------------------------------------------------------------
   TestGen() - Test finite field generator
   -------------------------------------------------------------------------- */

static void TestGen1()
{
   FEL a, b;
   int i;

   a = FfGen;
   b = a;
   for (i = 1; i < (int)FfOrder - 1; ++i) {
      if (b == FF_ONE) {
         Error("Generator has order %d",i);
      }
      b = FfMul(a,b);
   }
   if (b != FF_ONE) {
      Error("g^%d = %d, should be %d",b,FF_ONE);
   }
}


void TestGen(unsigned flags)
{
   while (NextField() > 0) {
      TestGen1();
   }
   flags = 0;
}


/* --------------------------------------------------------------------------
   TestInsertExtract() - Test vector insert/extract operations
   -------------------------------------------------------------------------- */

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
               Error("Read %d %d %d, expected %d %d %d",
                     FfExtract(x,pos + 1 - 1),FfExtract(x,pos + 2 - 1),
                     FfExtract(x,pos + 3 - 1),FTab[i1],FTab[i2],FTab[i3]);
            }

         }
      }
   }
}


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


void TestInsertExtract(unsigned flags)
{
   while (NextField() > 0) {
      TestInsertExtract1();
   }
   flags = 0;
}


/* --------------------------------------------------------------------------
   TestFindPiv() - Test FfFindPivot()
   -------------------------------------------------------------------------- */

static void TestFindPiv2(PTR x, int noc)
{
   FEL a;
   int pos;
   int i;

   if ((pos = FfFindPivot(x,&a)) >= 0) {
      Error("Found pivot %d at column %d in empty row",a,pos);
   }
   for (i = 0; i < noc; ++i) {
      FfInsert(x,i,FTab[i % (FfOrder - 1) + 1]);
   }
   for (i = 0; i < noc; ++i) {
      int k;
      for (k = 1; k < FfOrder; ++k) {
         FfInsert(x,i,FTab[k]);
         if (((pos = FfFindPivot(x,&a)) != i) || (a != FTab[k])) {
            Error("Found Pivot %d at %d, expected %d at %d",
                  a,pos,FTab[k],i);
         }
      }
      FfInsert(x,i,FF_ZERO);
   }
}


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
         Error("Found pivot at %d, expected %d",pos,i);
      }
      FfSetNoc(i);
      pos = FfFindPivot(x,&a);
      if (pos != -1) {
         Error("Found pivot at %d, expected %d",pos,-1);
      }
   }
}


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


void TestFindPiv(unsigned flags)
{
   while (NextField() > 0) {
      TestFindPiv1();
   }
   flags = 0;
}


/* --------------------------------------------------------------------------
   TestFelToInt() - Test FEL <--> integer conversion
   -------------------------------------------------------------------------- */

static void TestFelToInt1()
{
   long l, m;
   FEL a;

   if (FfFromInt(0) != FF_ZERO) {
      Error("FfFromInt(0) = %d, expected %d",FfFromInt(0),FF_ZERO);
   }
   if (FfFromInt(1) != FF_ONE) {
      Error("FfFromInt(1) = %d, expected %d",FfFromInt(1),FF_ONE);
   }
   if (FfToInt(FF_ZERO) != 0) {
      Error("FfToInt(%d) = %d, expected 0",FF_ZERO,0);
   }
   if (FfToInt(FF_ONE) != 1) {
      Error("FfToInt(%d) = %d, expected 1",FF_ONE,1);
   }
   for (l = 0; l < FfOrder; ++l) {
      FTab[l] = FfFromInt(l);
      if (!ISFEL(FTab[l])) {
         Error("FfFromInt(%d)=%d illegal",l,FTab[l]);
      }
      for (m = 0; m < l; ++m) {
         if (FTab[m] == FTab[l]) {
            Error("FfFromInt(%d) = zitof(%d)",m,l);
         }
      }
      if (FfToInt(FTab[l]) != l) {
         Error("FfToInt(FfFromInt(%d)=%d",l,FfToInt(FTab[l]));
      }
   }

   for (a = FF_ZERO, l = 0; l < FfChar; ++l, a = FfAdd(a,FF_ONE)) {
      if (FTab[l] != a) {
         Error("FfFromInt(%d)=%d, expected %d",FTab[l],a);
      }
   }

}


void TestFelToInt(unsigned flags)
{
   while (NextField() > 0) {
      TestFelToInt1();
   }
   flags = 0;
}


/* --------------------------------------------------------------------------
   TestSubfield() - Test subfield embedding/restriction
   -------------------------------------------------------------------------- */

static void TestSubfield1(int fld, int sub)
{
   FEL tabsub[256];
   FEL tabemb[256];
   int i;

   FfSetField(fld);
   FfSetField(sub);
   for (i = 0; i < sub; ++i) {
      tabsub[i] = FfFromInt(i);
   }
   FfSetField(fld);
   for (i = 0; i < sub; ++i) {
      tabemb[i] = FfEmbed(tabsub[i],sub);
      if (!ISFEL(tabemb[i])) {
         Error("Embedding of %d from GF(%d) to GF(%d) invalid",i,sub,fld);
      }
      if (FfRestrict(tabemb[i],sub) != tabsub[i]) {
         Error("FfRestrict(FfEmbed(%d))=%d",tabsub[i],
               FfRestrict(tabemb[i],sub));
      }
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
            Error("Embedding of GF(%d) into GF(%d) not closed (+)",fld,sub);
         }
         for (l = 0; l < sub && tabemb[l] != y; ++l) {
         }
         if (l >= sub) {
            Error("Embedding of GF(%d) into GF(%d) not closed (*)",fld,sub);
         }
      }
   }
}


void TestSubfields(unsigned flags)
{
   printf(" 2");
   fflush(stdout);
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

   printf(" 3");
   fflush(stdout);
   TestSubfield1(243,3);
   TestSubfield1(81,3);
   TestSubfield1(81,9);
   TestSubfield1(9,3);

   printf(" 5");
   fflush(stdout);
   TestSubfield1(125,5);
   TestSubfield1(25,5);

   printf(" 7");
   fflush(stdout);
   TestSubfield1(49,7);

   printf(" 11");
   fflush(stdout);
   TestSubfield1(121,11);
   flags = 0;
}


/* --------------------------------------------------------------------------
   TestRowOps() - Test row operations
   -------------------------------------------------------------------------- */

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
         Error("FfAddRow failed at %d+%d in column %d",
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
            Error("FfAddMulRow failed at %d+%d*%d in column %d",
                  FTab[(i + d1) % FfOrder],FTab[k],FTab[(i + d2) % FfOrder],i);
         }
      }
   }
}


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
               Error("FfAddrow failed at col %d: %d+%d=%d, expected %d",
                     i,FTab[ai],FTab[bi],ist,soll);
            }
         }
      }
   }
}


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


void TestRowOps(unsigned flags)
{
   while (NextField() > 0) {
      TestAddRow1();
   }
   flags = 0;
}


/* --------------------------------------------------------------------------
   TestMulRow() - Test FfMulRow()
   -------------------------------------------------------------------------- */

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
            Error("FfMulRow() failed");
         }
         if ((val /= FfOrder) == 0) {
            val = i;
         }
      }
   }

   FfMulRow(x,FF_ZERO);
   for (i = 0; i < noc; ++i) {
      if (FfExtract(x,i) != FF_ZERO) {
         Error("FfMulRow(x,0) failed");
      }
   }

}


void TestMulRow(unsigned flags)
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
   flags = 0;
}


/* --------------------------------------------------------------------------
   TestPtr() - Test row pointer arithmetic
   -------------------------------------------------------------------------- */

static void TestPtr2(PTR x)
{
   PTR xx;
   xx = x;
   FfStepPtr(&xx);
   if ((char *)xx != (char *)x + FfCurrentRowSize) {
      Error("Pointer increment error");
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


void TestPtr(unsigned flags)
{
   while (NextField() > 0) {
      TestPtr1();
   }
   flags = 0;
}


/* --------------------------------------------------------------------------
   TestRowSize() - Test row sizes
   -------------------------------------------------------------------------- */

static void TestRowSize1()
{
   int noc;
   for (noc = 0; noc < 100; ++noc) {
      int rs;
      int diff;
      rs = FfRowSize(noc);
      if ((rs < 0) || (rs > noc + sizeof(long))) {
         Error("FfRowSize(%d) = %d out of range",noc,rs);
      }
      diff = rs - FfTrueRowSize(noc);
      if ((diff < 0) || (diff >= sizeof(long))) {
         Error("FfRowSize() and FfTrueRowSize() differ too much");
      }
   }
}


void TestRowSize(unsigned flags)
{
   while (NextField() > 0) {
      TestRowSize1();
   }
   flags = 0;
}


/* --------------------------------------------------------------------------
   TestCmpRows() - Test FfCmpRows()
   -------------------------------------------------------------------------- */

static void TestCmpRows2(PTR m1, PTR m2, int noc)
{
   int i;

   FfMulRow(m1,FF_ZERO);
   FfMulRow(m2,FF_ZERO);

   for (i = 1; i < FfOrder; ++i) {
      int k;
      for (k = 0; k < noc; ++k) {
         if (FfCmpRows(m2,m1) != 0) {
            Error("Rows are different");
         }
         FfInsert(m1,k,FTab[i]);
         if (FfCmpRows(m2,m1) == 0) {
            Error("Rows are still equal");
         }
         FfInsert(m2,k,FTab[i]);
         if (FfCmpRows(m2,m1) != 0) {
            Error("Rows are still different");
         }
      }
   }
}


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


void TestCmpRows(unsigned flags)
{
   while (NextField() > 0) {
      TestCmpRows1();
   }
   flags = 0;
}