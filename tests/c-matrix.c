////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check functions for matrices.
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

Matrix_t *RndMat(int fl, int nor, int noc)
{
   Matrix_t *a = MatAlloc(fl,nor,noc);
   int i;
   for (i = 0; i < nor; ++i) {
      PTR ax = MatGetPtr(a,i);
      int k;
      for (k = 0; k < noc; ++k) {
         FfInsert(ax,k,FTab[MtxRandomInt(FfOrder)]);
      }
   }
   return a;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

#define NMAT 5

static void TestMatAlloc1(int fl)
{
   static const int nor[NMAT] = { 0,0,1,1,9 };
   static const int noc[NMAT] = { 0,1,0,1,9 };
   Matrix_t *m[NMAT];
   int i;

   for (i = 0; i < NMAT; ++i) {
      m[i] = MatAlloc(fl,nor[i],noc[i]);
   }
   for (i = 0; i < NMAT; ++i) {
      MatIsValid(m[i]);
      ASSERT(m[i]->Field == fl);
      ASSERT(m[i]->Nor == nor[i]);
      ASSERT(m[i]->Noc == noc[i]);
   }
   for (i = 0; i < NMAT; ++i) {
      ASSERT_EQ_INT(MatFree(m[i]), 0);
   }
   TstStartErrorChecking();
   for (i = 0; i < NMAT; ++i) {
      ASSERT(!MatIsValid(m[i]));
      ASSERT(TstHasError());
   }
   TstStopErrorChecking();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F MatrixAllocation(unsigned flags)
{
   while (NextField() > 0) {
      TestMatAlloc1(FfOrder);
   }
   flags = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void ChkEch(Matrix_t *mat)
{
   int i;

   for (i = 0; i < mat->Nor; ++i) {
      PTR p = MatGetPtr(mat,i);
      int k;
      FEL f;

      ASSERT_EQ_INT(FfFindPivot(p,&f), mat->PivotTable[i]);
      for (k = 0; k < i; ++k) {
	 ASSERT_EQ_INT(FfExtract(p,mat->PivotTable[k]), FF_ZERO);
      }
   }
   for (i = mat->Nor; i < mat->Noc; ++i) {
      int k;
      int piv = mat->PivotTable[i];
      ASSERT(piv >= 0);
      ASSERT(piv < mat->Noc);
      for (k = 0; k < i; ++k) {
	 ASSERT(mat->PivotTable[k] != piv);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMatEchelonize1(Matrix_t *m, int size)
{
   int i;

   for (i = 0; i < size; ++i) {
      PTR p = MatGetPtr(m,i);
      int k;
      FfMulRow(p,FF_ZERO);
      for (k = size - i - 1; k < size; ++k) {
         FfInsert(p,k,FF_ONE);
      }
   }
   ASSERT_EQ_INT(MatEchelonize(m), size);
   ChkEch(m);
   for (i = 0; i < size; ++i) {
      PTR p = MatGetPtr(m,i);
      int k;
      ASSERT_EQ_INT(m->PivotTable[i], size - i - 1);
      for (k = 0; k < size; ++k) {
         FEL f = FfExtract(p,k);
         if ((f == FF_ZERO) ^ (k != size - i - 1)) {
            TST_FAIL2("Unexpected element %d at column %d",f,k);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMatEchelonize2(Matrix_t *m, int size)
{
   int i;
   static unsigned long x = 0;

   for (i = 1; i <= size; ++i) {
      int k;
      PTR p = MatGetPtr(m,i - 1);
      for (k = 0; k < size; ++k) {
         FfInsert(p,k,FTab[(x >> 3) % FfOrder]);
         x = 69069 * x + 3;
      }
   }
   ASSERT(MatEchelonize(m) >= 5);
   ChkEch(m);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F TestMatEchelonize()
{
   const int size = 10;

   while (NextField() > 0) {
      Matrix_t *m = MatAlloc(FfOrder,size,size);
      TestMatEchelonize1(m,size);
      TestMatEchelonize2(m,size);
      MatFree(m);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMatCompare1(Matrix_t *a, Matrix_t *b, int size)
{
   ASSERT(MatCompare(a,b) == 0);
   ASSERT(MatCompare(b,a) == 0);
   for (int i = 0; i < size; ++i) {
      PTR pa = MatGetPtr(a,i);
      PTR pb = MatGetPtr(b,i);
      FfInsert(pa,0,FF_ONE);
      ASSERT(MatCompare(a,b) > 0);
      ASSERT(MatCompare(b,a) < 0);

      FfInsert(pb,0,FF_ONE);
      ASSERT(MatCompare(a,b) == 0);
      ASSERT(MatCompare(b,a) == 0);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void Check2(int fld1, int nor1, int noc1, int fld2, int nor2,
                   int noc2)
{
   Matrix_t *a, *b;
   a = MatAlloc(fld1,nor1,noc1);
   b = MatAlloc(fld2,nor2,noc2);
   ASSERT(MatCompare(a,b) < 0);
   ASSERT(MatCompare(b,a) > 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMatCompare2()
{
   Check2(2,20,20,3,10,10);
   Check2(2,20,20,2,10,30);
   Check2(2,20,20,2,30,20);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F MatrixCompare()
{

   while (NextField() > 0) {
      int size;
      for (size = 2; size < 10; ++size) {
         Matrix_t *a = MatAlloc(FfOrder,size,size);
         Matrix_t *b = MatAlloc(FfOrder,size,size);
         TestMatCompare1(a,b,size);
         MatFree(a);
         MatFree(b);
      }
   }
   TestMatCompare2();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMatClean1()
{
   Matrix_t *a =
      MkMat(4,6, 1,0,0,0,0,0,  0,1,1,0,0,0, 0,0,0,0,1,0, 0,0,1,1,0,0);
   Matrix_t *b =
      MkMat(4,6, 0,0,0,0,0,1,  0,0,0,0,0,0, 1,1,1,1,1,1, 1,0,1,0,1,0);
   Matrix_t *c =
      MkMat(2,6, 0,0,0,0,0,1,  0,0,0,1,0,0);

   MatEchelonize(a);
   if (MatClean(b,a) != 2) {
      TST_FAIL("MatClean() failed");
   }
   if (MatCompare(b,c) != 0) {
      TST_FAIL("MatClean() produced wrong result");
   }

   MatFree(a);
   MatFree(b);
   MatFree(c);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F MatrixClean()
{
   while (NextField() > 0) {
      TestMatClean1();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMatInv1()
{
   int i;

   for (i = 0; i < 20; ++i) {
      Matrix_t *a = MatId(FfOrder,i);
      Matrix_t *ai = MatInverse(a);
      if (MatCompare(a,ai) != 0) {
         TST_FAIL("Wrong inverse of identity matrix");
      }
      MatFree(a);
      MatFree(ai);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMatInv2()
{
   Matrix_t *a =
      MkMat(5,5, 1,2,3,0,2,  0,0,0,1,1, 0,0,1,1,0, 0,1,2,3,0, 0,0,0,0,1);
   Matrix_t *ai, *id;
   ai = MatInverse(a);
   MatMul(ai,a);
   id = MatId(FfOrder,5);
   if (MatCompare(ai,id) != 0) {
      TST_FAIL("M * 1/M != 1");
   }
   MatFree(a);
   MatFree(ai);
   MatFree(id);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F TestMatrixInversion()
{
   while (NextField() > 0) {
      TestMatInv1();
      TestMatInv2();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMatDup1(int fl, int nor, int noc)
{
   Matrix_t *a, *b;

   a = RndMat(fl,nor,noc);
   b = MatDup(a);
   if (MatCompare(a,b) != 0) {
      TST_FAIL("MatDup() failed");
   }
   MatFree(a);
   MatFree(b);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F MatrixDuplication()
{
   int nor, noc;
   while (NextField() > 0) {
      for (nor = 0; nor < 10; ++nor) {
         for (noc = 0; noc < 10; ++noc) {
            TestMatDup1(FfOrder,nor,noc);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestNullSpace1(int fl, int dim)
{
   Matrix_t *a, *b;
   a = MatId(fl,dim);
   b = MatNullSpace(a);
   if (!MatIsValid(b) || (b->Noc != dim) || (b->Nor != 0)) {
      TST_FAIL("NullSpace(Id) failed");
   }
   MatFree(b);
   b = MatNullSpace__(MatAlloc(fl,dim,dim));
   if (MatCompare(a,b) != 0) {
      TST_FAIL("NullSpace(0) != Id");
   }
   MatFree(a);
   MatFree(b);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestNullSpace2(int fl, int dim)
{
   Matrix_t *a, *b;
   int i;

   a = RndMat(fl,dim + 3,dim);
   b = MatNullSpace(a);
   if (b->Nor < 3) {
      TST_FAIL("Unexpected null-space dimension");
   }
   MatMul(b,a);
   FfSetNoc(b->Noc);
   for (i = 0; i < b->Nor; ++i) {
      FEL f;
      if (FfFindPivot(MatGetPtr(b,i),&f) >= 0) {
         TST_FAIL("Wrong null-space");
      }
   }
   MatFree(a);
   MatFree(b);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F NullSpace()
{
   int nor;

   while (NextField() > 0) {
      for (nor = 0; nor < 10; ++nor) {
         TestNullSpace1(FfOrder,nor);
         if (nor > 0) {
            TestNullSpace2(FfOrder,nor);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMatOrder2(Matrix_t *a, int order)
{
   ASSERT_EQ_INT(MatOrder(a), order);
   MatFree(a);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMatOrder1()
{
   Matrix_t *a;

   a = MkMat(5,5, 1,0,0,0,0,  0,1,0,0,0, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,1);
   TestMatOrder2(a,1);

   a = MkMat(3,3, -1,1,0, -1,0,1, 0,0,1);
   TestMatOrder2(a,3);

   a = MkMat(5,5, 0,1,0,-1,0,  1,1,0,-1,1, -1,1,0,0,0, 0,1,0,-1,1, -1,0,1,0,0);
   TestMatOrder2(a,6);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F MatrixOrder()
{
   while (NextField() > 0) {
      TestMatOrder1();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMatCut1(int fl)
{
   Matrix_t *a, *b;
   int anor = 10, anoc = 20;
   int i;

   a = RndMat(fl,anor,anoc);

   b = MatCut(a,0,0,-1,-1);
   if (MatCompare(a,b) != 0) {
      TST_FAIL("MatCut(...0,0,-1,-1) failed");
   }
   MatFree(b);
   b = MatCutRows(a,0,-1);
   if (MatCompare(a,b) != 0) {
      TST_FAIL("MatCutRows(...,0,-1) failed");
   }
   MatFree(b);

   for (i = 0; i < anor * anoc * 10; ++i) {
      int bnor, bnoc;
      int nor0, noc0;
      int r;

      nor0 = MtxRandomInt(anor);
      noc0 = MtxRandomInt(anoc);
      bnor = MtxRandomInt(anor - nor0);
      bnoc = MtxRandomInt(anoc - noc0);
      b = MatCut(a,nor0,noc0,bnor,bnoc);
      for (r = 0; r < bnor; ++r) {
         int c;
         for (c = 0; c < bnoc; ++c) {
            FEL fa, fb;
            FfSetNoc(anoc);
            fa = FfExtract(MatGetPtr(a,nor0 + r),noc0 + c);
            FfSetNoc(bnoc);
            fb = FfExtract(MatGetPtr(b,r),c);
            if (fa != fb) {
               TST_FAIL6("Error at (%d+%d,%d+%d): %d != %d",nor0,r,noc0,c,
                     (int)fa,(int)fb);
            }
         }
      }
      MatFree(b);
   }
   MatFree(a);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F TestMatrixCut()
{
   while (NextField() > 0) {
      TestMatCut1(FfOrder);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMatCopy1(int fl)
{
   Matrix_t *a, *b;
   int nor = 10, noc = 20;
   int i;

   a = RndMat(fl,nor,noc);
   b = MatAlloc(fl,nor,noc);

   for (i = 0; i < nor * noc * 10; ++i) {
      int sr0, sc0, snor, snoc;
      int dr0, dc0;
      int r;

      sr0 = MtxRandomInt(nor);
      sc0 = MtxRandomInt(noc);
      snor = MtxRandomInt(nor - sr0);
      snoc = MtxRandomInt(noc - sc0);
      dr0 = MtxRandomInt(nor - snor);
      dc0 = MtxRandomInt(noc - snoc);

      MatMulScalar(b,FF_ZERO);
      if (MatCopyRegion(b,dr0,dc0,a,sr0,sc0,snor,snoc) != 0) {
         TST_FAIL("MatCopyRegion() failed");
      }
      for (r = 0; r < nor; ++r) {
         int c;
         for (c = 0; c < noc; ++c) {
            FEL fb = FfExtract(MatGetPtr(b,r),c);
            if ((r < dr0) || (r >= dr0 + snor) || (c < dc0) || (c >= dc0 + snoc)) {
               if (fb != FF_ZERO) {
                  TST_FAIL3("Found %d at (%d,%d), expected 0",fb,r,c);
               }
            } else {
	       ASSERT_EQ_INT(FfExtract(MatGetPtr(a,sr0 + r - dr0),sc0 + c - dc0), fb);
            }
         }
      }
   }
   MatFree(a);
   MatFree(b);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F MatrixCopy(unsigned flags)
{
   while (NextField() > 0) {
      TestMatCopy1(FfOrder);
   }
   flags = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMatAddMul2(int fl, int nor, int noc, Matrix_t *a, Matrix_t *b,
                           Matrix_t *c)
{
   int i;
   for (i = 0; i < fl; i += i / 10 + 1) {
      int r;
      FEL f = FfFromInt(i);
      MatCopyRegion(b,0,0,a,0,0,-1,-1);
      MatAddMul(a,c,f);
      for (r = 0; r < nor; ++r) {
         int co;
         for (co = 0; co < noc; ++co) {
	    FEL fa = FfExtract(MatGetPtr(a,r),co);
            FEL fb = FfExtract(MatGetPtr(b,r),co);
            FEL fc = FfExtract(MatGetPtr(c,r),co);
            fb = FfAdd(fb,FfMul(fc,f));
	    ASSERT_EQ_INT(fa, fb);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMatAddMul1(int fl)
{
   int nor;
   for (nor = 0; nor < 20; nor += nor / 5 + 1) {
      int noc;
      for (noc = 0; noc < 20; noc += noc / 5 + 1) {
         Matrix_t *a = RndMat(fl,nor,noc);
         Matrix_t *b = MatDup(a);
         Matrix_t *c = RndMat(fl,nor,noc);
         TestMatAddMul2(fl,nor,noc,a,b,c);
         MatFree(a);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F MatrixMultiplyAdd()
{
   while (NextField() > 0) {
      TestMatAddMul1(FfOrder);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMatId2(int fl, int dim)
{
   Matrix_t *m;
   int i;

   m = MatId(fl,dim);
   ASSERT_EQ_INT(m->Field, fl);
   ASSERT_EQ_INT(m->Nor, dim);
   ASSERT_EQ_INT(m->Noc, dim);
   for (i = 0; i < dim; ++i) {
      int k;
      PTR r = MatGetPtr(m,i);
      for (k = 0; k < dim; ++k) {
         const FEL f = FfExtract(r,k);
         if (k == i) ASSERT_EQ_INT(f, FF_ONE); else ASSERT_EQ_INT(f, FF_ZERO);
      }
   }
   MatFree(m);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F MatrixIdentity()
{
   while (NextField() > 0) {
       for (int dim = 0; dim < 100; dim += dim / 10 + 1) {
	   TestMatId2(FfOrder,dim);
       }
   }
}
