////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check functions for matrices.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

Matrix_t *RndMat(int fl, int nor, int noc)
{
   Matrix_t *a = matAlloc(fl,nor,noc);
   int i;
   for (i = 0; i < nor; ++i) {
      PTR ax = matGetPtr(a,i);
      int k;
      for (k = 0; k < noc; ++k) {
         ffInsert(ax,k,FTab[mtxRandomInt(ffOrder)]);
      }
   }
   return a;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

#define NMAT 5

TstResult Matrix_Allocation(int q)
{
   static const int nor[NMAT] = { 0,0,1,1,9 };
   static const int noc[NMAT] = { 0,1,0,1,9 };
   Matrix_t *m[NMAT];
   int i;

   for (i = 0; i < NMAT; ++i) {
      m[i] = matAlloc(ffOrder,nor[i],noc[i]);
   }
   for (i = 0; i < NMAT; ++i) {
      ASSERT(matIsValid(m[i]));
      ASSERT(m[i]->Field == ffOrder);
      ASSERT(m[i]->Nor == nor[i]);
      ASSERT(m[i]->Noc == noc[i]);
   }
   for (i = 0; i < NMAT; ++i) {
      ASSERT_EQ_INT(matFree(m[i]), 0);
      ASSERT(!matIsValid(m[i]));
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Matrix_ThrowsOnDoubleFree()
{
   Matrix_t *m = matAlloc(3, 20, 30);
   ASSERT(matIsValid(m));
   matFree(m);
   ASSERT(!matIsValid(m));
   ASSERT_ABORT(matFree(m));
   return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////

static int ChkEch(Matrix_t *mat)
{
   for (uint32_t i = 0; i < mat->Nor; ++i) {
      PTR p = matGetPtr(mat,i);
      FEL f;

      ASSERT_EQ_INT(ffFindPivot(p,&f,mat->Noc), mat->PivotTable[i]);
      for (uint32_t k = 0; k < i; ++k) {
	 ASSERT_EQ_INT(ffExtract(p,mat->PivotTable[k]), FF_ZERO);
      }
   }
   for (uint32_t i = mat->Nor; i < mat->Noc; ++i) {
      uint32_t piv = mat->PivotTable[i];
      ASSERT(piv < mat->Noc);
      for (uint32_t k = 0; k < i; ++k) {
	 ASSERT(mat->PivotTable[k] != piv);
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestMatEchelonize1(Matrix_t *m, int size)
{
   int i;

   for (i = 0; i < size; ++i) {
      PTR p = matGetPtr(m,i);
      int k;
      ffMulRow(p,FF_ZERO, size);
      for (k = size - i - 1; k < size; ++k) {
         ffInsert(p,k,FF_ONE);
      }
   }
   ASSERT_EQ_INT(matEchelonize(m), size);
   ASSERT(ChkEch(m) == 0);
   for (i = 0; i < size; ++i) {
      PTR p = matGetPtr(m,i);
      int k;
      ASSERT_EQ_INT(m->PivotTable[i], size - i - 1);
      for (k = 0; k < size; ++k) {
         FEL f = ffExtract(p,k);
         ASSERT(!((f == FF_ZERO) ^ (k != size - i - 1)));
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestMatEchelonize2(Matrix_t *m, int size)
{
   int i;
   static unsigned long x = 0;

   for (i = 1; i <= size; ++i) {
      int k;
      PTR p = matGetPtr(m,i - 1);
      for (k = 0; k < size; ++k) {
         ffInsert(p,k,FTab[(x >> 3) % ffOrder]);
         x = 69069 * x + 3;
      }
   }
   int result = 0;
   ASSERT(matEchelonize(m) >= 5);
   result |= ChkEch(m);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Matrix_Echelonize(int q)
{
   const int size = 10;
   Matrix_t *m = matAlloc(ffOrder,size,size);
   int result = 0;
   result |= TestMatEchelonize1(m,size);
   result |= TestMatEchelonize2(m,size);
   matFree(m);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestMatCompare1(Matrix_t *a, Matrix_t *b, int size)
{
   ASSERT(matCompare(a,b) == 0);
   ASSERT(matCompare(b,a) == 0);
   for (int i = 0; i < size; ++i) {
      PTR pa = matGetPtr(a,i);
      PTR pb = matGetPtr(b,i);
      ffInsert(pa,0,FF_ONE);
      ASSERT(matCompare(a,b) != 0);
      ASSERT(matCompare(b,a) != 0);

      ffInsert(pb,0,FF_ONE);
      ASSERT(matCompare(a,b) == 0);
      ASSERT(matCompare(b,a) == 0);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Matrix_Compare1(int q)
{
   int result = 0;
   for (int size = 2; size < 10; ++size) {
      Matrix_t *a = matAlloc(ffOrder,size,size);
      Matrix_t *b = matAlloc(ffOrder,size,size);
      result |= TestMatCompare1(a,b,size);
      matFree(a);
      matFree(b);
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int check2(int nor1, int noc1, int nor2, int noc2, int expectedResult)
{
   Matrix_t *a, *b;
   a = matAlloc(ffOrder,nor1,noc1);
   b = matAlloc(ffOrder,nor2,noc2);
   ASSERT_EQ_INT(matCompare(a,b), expectedResult);
   ASSERT_EQ_INT(matCompare(b,a), -expectedResult);
   matFree(a);
   matFree(b);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Matrix_CompareSize(int q)
{
   int result = 0;
   for (int n = 1; result == 0 && n < 16; ++n) {
      result |= check2(n, n, n, n+1, -1);
      result |= check2(n, n, n, n-1, 1);
      result |= check2(n, n, n+1, n - 1, 1);
      result |= check2(n, n, n-1, n + 1, -1);
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Matrix_Clean(int q)
{
   Matrix_t *a = MkMat(4,6, 1,0,0,0,0,0,  0,1,1,0,0,0, 0,0,0,0,1,0, 0,0,1,1,0,0);
   Matrix_t *b = MkMat(4,6, 0,0,0,0,0,1,  0,0,0,0,0,0, 1,1,1,1,1,1, 1,0,1,0,1,0);
   Matrix_t *c = MkMat(2,6, 0,0,0,0,0,1,  0,0,0,1,0,0);

   matEchelonize(a);
   ASSERT(matClean(b,a) == 2);
   ASSERT(matCompare(b,c) == 0);

   matFree(a);
   matFree(b);
   matFree(c);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestMatInv1()
{
   for (int i = 0; i < 20; ++i) {
      Matrix_t *a = matId(ffOrder,i);
      Matrix_t *ai = matInverse(a);
      ASSERT(matCompare(a,ai) == 0);
      matFree(a);
      matFree(ai);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestMatInv2()
{
   Matrix_t *a = MkMat(5,5, 1,2,3,0,2,  0,0,0,1,1, 0,0,1,1,0, 0,1,2,3,0, 0,0,0,0,1);
   Matrix_t *ai, *id;
   ai = matInverse(a);
   matMul(ai,a);
   id = matId(ffOrder,5);
   ASSERT(matCompare(ai,id) == 0);
   matFree(a);
   matFree(ai);
   matFree(id);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Matrix_Inversion(int q)
{
   int result = 0;
   result |= TestMatInv1();
   result |= TestMatInv2();
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Matrix_Duplication(int q)
{
   for (int nor = 0; nor < 10; ++nor) {
      for (int noc = 0; noc < 10; ++noc) {
         Matrix_t* a = RndMat(ffOrder,nor,noc);
         Matrix_t* b = matDup(a);
         ASSERT(matCompare(a,b) == 0);
         matFree(a);
         matFree(b);
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestNullSpace1(int dim)
{
   Matrix_t* a = matId(ffOrder,dim);
   Matrix_t* b = matNullSpace(a);
   ASSERT(matIsValid(b));
   ASSERT(b->Noc == dim);
   ASSERT(b->Nor == 0);
   matFree(b);

   b = matNullSpace__(matAlloc(ffOrder,dim,dim));
   ASSERT(matCompare(a,b) == 0);
   matFree(a);
   matFree(b);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestNullSpace2(int dim)
{
   Matrix_t* a = RndMat(ffOrder,dim + 3,dim);
   Matrix_t* b = matNullSpace(a);
   ASSERT(b->Nor >= 3);
   matMul(b,a);
   for (int i = 0; i < b->Nor; ++i) {
      FEL f;
      ASSERT(ffFindPivot(matGetPtr(b,i),&f, b->Noc) == MTX_NVAL);
   }
   matFree(a);
   matFree(b);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Matrix_NullSpace(int q)
{
   int result = 0;
   for (int nor = 0; result == 0 && nor < 10; ++nor) {
      result |= TestNullSpace1(nor);
      if (nor > 0) {
         result |= TestNullSpace2(nor);
      }
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Matrix_Order(int q)
{
   Matrix_t *a;

   a = MkMat(5,5, 1,0,0,0,0,  0,1,0,0,0, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,1);
   ASSERT_EQ_INT(matOrder(a), 1);
   matFree(a);

   a = MkMat(3,3, -1,1,0, -1,0,1, 0,0,1);
   ASSERT_EQ_INT(matOrder(a), 3);
   matFree(a);

   a = MkMat(5,5, 0,1,0,-1,0,  1,1,0,-1,1, -1,1,0,0,0, 0,1,0,-1,1, -1,0,1,0,0);
   ASSERT_EQ_INT(matOrder(a), 6);
   matFree(a);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Matrix_Cut(int q)
{
   Matrix_t *a, *b;
   int anor = 10, anoc = 20;
   int i;

   a = RndMat(ffOrder,anor,anoc);

   b = matCut(a,0,0,-1,-1);
   ASSERT(matCompare(a,b) == 0);
   matFree(b);
   b = matCutRows(a,0,-1);
   ASSERT(matCompare(a,b) == 0);
   matFree(b);

   for (i = 0; i < anor * anoc * 10; ++i) {
      int bnor, bnoc;
      int nor0, noc0;
      int r;

      nor0 = mtxRandomInt(anor);
      noc0 = mtxRandomInt(anoc);
      bnor = mtxRandomInt(anor - nor0);
      bnoc = mtxRandomInt(anoc - noc0);
      b = matCut(a,nor0,noc0,bnor,bnoc);
      for (r = 0; r < bnor; ++r) {
         int c;
         for (c = 0; c < bnoc; ++c) {
            FEL fa, fb;
            fa = ffExtract(matGetPtr(a,nor0 + r),noc0 + c);
            fb = ffExtract(matGetPtr(b,r),c);
            ASSERT(fa == fb);
         }
      }
      matFree(b);
   }
   matFree(a);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Matrix_Copy(int q)
{
   Matrix_t *a, *b;
   int nor = 10, noc = 20;
   int i;

   a = RndMat(ffOrder,nor,noc);
   b = matAlloc(ffOrder,nor,noc);

   for (i = 0; i < nor * noc * 10; ++i) {
      int sr0 = mtxRandomInt(nor);
      int sc0 = mtxRandomInt(noc);
      int snor = mtxRandomInt(nor - sr0);
      int snoc = mtxRandomInt(noc - sc0);
      int dr0 = mtxRandomInt(nor - snor);
      int dc0 = mtxRandomInt(noc - snoc);

      matMulScalar(b,FF_ZERO);
      ASSERT(matCopyRegion(b,dr0,dc0,a,sr0,sc0,snor,snoc) == 0);
      for (int r = 0; r < nor; ++r) {
         for (int c = 0; c < noc; ++c) {
            FEL fb = ffExtract(matGetPtr(b,r),c);
            if ((r < dr0) || (r >= dr0 + snor) || (c < dc0) || (c >= dc0 + snoc)) {
		ASSERT_EQ_INT(fb, FF_ZERO);
            } else {
	       ASSERT_EQ_INT(ffExtract(matGetPtr(a,sr0 + r - dr0),sc0 + c - dc0), fb);
            }
         }
      }
   }
   matFree(a);
   matFree(b);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestMatAddMul2(int nor, int noc, Matrix_t *a, Matrix_t *b, Matrix_t *c)
{
   for (int i = 0; i < ffOrder; i += i / 10 + 1) {
      FEL f = ffFromInt(i);
      matCopyRegion(b,0,0,a,0,0,-1,-1);
      matAddMul(a,c,f);
      for (int r = 0; r < nor; ++r) {
         int co;
         for (co = 0; co < noc; ++co) {
	    FEL fa = ffExtract(matGetPtr(a,r),co);
            FEL fb = ffExtract(matGetPtr(b,r),co);
            FEL fc = ffExtract(matGetPtr(c,r),co);
            fb = ffAdd(fb,ffMul(fc,f));
	    ASSERT_EQ_INT(fa, fb);
         }
      }
   }
   return 0;
}

TstResult Matrix_MultiplyAdd(int q)
{
   int result = 0;
   for (int nor = 0; nor < 20; nor += nor / 5 + 1) {
      for (int noc = 0; noc < 20; noc += noc / 5 + 1) {
         Matrix_t *a = RndMat(ffOrder,nor,noc);
         Matrix_t *b = matDup(a);
         Matrix_t *c = RndMat(ffOrder,nor,noc);
         result |= TestMatAddMul2(nor,noc,a,b,c);
         matFree(a);
      }
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestMatId2(int fl, int dim)
{
   Matrix_t *m;
   int i;

   m = matId(fl,dim);
   ASSERT_EQ_INT(m->Field, fl);
   ASSERT_EQ_INT(m->Nor, dim);
   ASSERT_EQ_INT(m->Noc, dim);
   for (i = 0; i < dim; ++i) {
      int k;
      PTR r = matGetPtr(m,i);
      for (k = 0; k < dim; ++k) {
         const FEL f = ffExtract(r,k);
         if (k == i)
	     ASSERT_EQ_INT(f, FF_ONE);
	 else 
	     ASSERT_EQ_INT(f, FF_ZERO);
      }
   }
   matFree(m);
   return 0;
}


TstResult Matrix_Identity(int q)
{
    int result = 0;
    for (int dim = 0; result == 0 && dim < 20; ++dim) {
	result |= TestMatId2(ffOrder,dim);
    }
    return result;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
