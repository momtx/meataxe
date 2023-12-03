////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for matrix functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestMapRow1(PTR mat, PTR a, PTR b, int noc)
{
   int i;

   /* Make <mat> the unit matrix.
      --------------------------- */
   for (i = 0; i < noc; ++i) {
      PTR m = (PTR)((char *)mat + ffSize(i, noc));
      ffInsert(m,i,FF_ONE);
   }

   for (i = 0; i < noc; ++i) {
      int k;

      // Check that the i-th basis vector is mapped to itself.
      ffMulRow(a,FF_ZERO, noc);
      ffInsert(a,i,FF_ONE);
      ffMapRow(b, a,mat,noc,noc);
      for (k = 0; k < noc; ++k) {
         ASSERT((ffExtract(b,k) == FF_ZERO) ^ (k == i));
      }

      
      ffMapRow(b, a,mat,i,noc);
      for (k = 0; k < noc; ++k) {
         ASSERT_EQ_INT(ffExtract(b,k), FF_ZERO);
      }
   }

   for (i = 0; i < noc; ++i) {
      ffInsert(a,i,FTab[i % ffOrder]);
   }
   ffMapRow(b, a,mat,noc,noc);
   for (i = 0; i < noc; ++i) {
      ASSERT_EQ_INT(ffExtract(b,i), FTab[i % ffOrder]);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_MapRow(int q)
{
   int size = 10;

   PTR mat, a, b;
   mat = ffAlloc(size, size);
   a = ffAlloc(1, size);
   b = ffAlloc(1, size);
   int result = TestMapRow1(mat,a,b,size);
   sysFree(mat);
   sysFree(a);
   sysFree(b);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestSumInter1(int noc)
{
   PTR w1, w2;
   uint32_t *piv;
   uint32_t nor1, nor2;
   uint32_t i, k;

   nor1 = nor2 = noc;
   w1 = ffAlloc(nor1 + nor2, noc);
   w2 = ffAlloc(nor1 + nor2, noc);
   piv = NALLOC(uint32_t,nor1 + nor2);
   k = 0;
   for (i = 0; i < nor1; ++i) {
      PTR x = ffGetPtr(w1,i,noc);
      ffInsert(x,k,FF_ONE);
      if (k % 3 == 0) {
         ++k;
      } else {
         k += 2;
      }
      if (k >= noc) {
         k = 0;
      }
   }
   k = 0;
   for (i = nor1; i < nor1 + nor2; ++i) {
      PTR x = ffGetPtr(w1,i,noc);
      ffInsert(x,k,FF_ONE);
      if (k % 3 == 0) {
         k += 2;
      } else {
         ++k;
      }
      if (k >= noc) {
         k = 0;
      }
   }

   //printrows("v", w1,nor1);
   //printrows("w", ffGetPtr(w1,nor1),nor2);

   ffSumAndIntersection(noc, w1,&nor1,&nor2,w2,piv);

   //printrows("sum", w1,nor1);
   //printrows("intersection", ffGetPtr(w2,nor1), nor2);

   ASSERT_EQ_INT(nor1, noc);	// whole space
   ASSERT_EQ_INT(nor2, (noc - 1) / 3 + 1);
   for (i = 0; i < nor2; ++i) {
      ASSERT_EQ_INT(piv[nor1 + i] % 3, 0);
   }

   sysFree(w1);
   sysFree(w2);
   sysFree(piv);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Checks that u is a subspace of v. As a side effect,
// - v is echelonized 
// - u is cleared (all rows set to zero).
static int CheckIsSubspace(uint32_t noc, PTR u, uint32_t udim, PTR v, uint32_t vdim)
{
   uint32_t *piv = NALLOC(uint32_t,vdim);
   uint32_t i, vrank;
   PTR x, y, row;

   //printf("CheckIsSubspace\n");
   //printrows("u",u,udim);
   //printrows("v",v,vdim);


   /* Echelonize v.
      ------------- */
   y = v;
   vrank = 0;
   for (i = 0, x = v; i < vdim; ++i, ffStepPtr(&x, noc)) {
      FEL f;
      ffCleanRow(x,v,vrank,noc,piv);
      uint32_t p = ffFindPivot(x,&f,noc);
      if (p == MTX_NVAL) {
         continue;
      }
      if (i > vrank) {
         ffCopyRow(y,x, noc);
      }
      piv[vrank++] = p;
      ffStepPtr(&y,noc);
   }
   //printrows("v chelonized", v, vrank);
   //printf("rank(v)=%d, cleaning u\n", vrank);

   /* Clean u with v.
      --------------- */
   row = ffAlloc(1, noc);
   for (i = 0, x = u; i < udim; ++i, ffStepPtr(&x,noc)) {
      FEL f;
      ffCopyRow(row,x, noc);
      ffCleanRow(row,v,vrank,noc,piv);
      ASSERT(ffFindPivot(row,&f, noc) == MTX_NVAL);
   }
   sysFree(row);
   sysFree(piv);
   //printf("CheckIsSubspace ok\n");
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestSumInter2(int noc)
{
   PTR v, w, wrk1, wrk2, x;
   uint32_t nor1, nor2;
   uint32_t vdim, wdim;
   uint32_t i, k;

   /* Allocate buffers.
      ----------------- */
   vdim = nor1 = mtxRandomInt(noc + 1);
   wdim = nor2 = mtxRandomInt(noc + 1);
   v = ffAlloc(vdim, noc);
   w = ffAlloc(wdim, noc);
   wrk1 = ffAlloc(nor1 + nor2, noc);
   wrk2 = ffAlloc(nor1 + nor2, noc);
   uint32_t *piv = NALLOC(uint32_t,nor1 + nor2);

   // Fill with random values.
   for (i = 0, x = v; i < nor1; ++i, ffStepPtr(&x, noc)) {
      for (k = 0; k < noc; ++k) {
         ffInsert(x,k,FTab[mtxRandomInt(ffOrder)]);
      }
   }
   for (i = 0, x = w; i < nor2; ++i, ffStepPtr(&x, noc)) {
      for (k = 0; k < noc; ++k) {
         ffInsert(x,k,FTab[mtxRandomInt(ffOrder)]);
      }
   }
   //printrows("v",v,vdim);
   //printrows("w",w,wdim);

   memcpy(wrk1,v,ffSize(nor1, noc));
   memcpy(ffGetPtr(wrk1,nor1,noc),w,ffSize(nor2,noc));
   ffSumAndIntersection(noc,wrk1,&nor1,&nor2,wrk2,piv);

   //printrows("sum", wrk1, nor1);
   //printrows("intersection", ffGetPtr(wrk2, nor1), nor2);
  

   /* Check relations between V, W, V+w, and VW.
      ------------------------------------------ */
   int result = 0;
   result |= CheckIsSubspace(noc, v,vdim,wrk1,nor1);
   result |= CheckIsSubspace(noc, w,wdim,wrk1,nor1);
   result |= CheckIsSubspace(noc, ffGetPtr(wrk2,nor1,noc),nor2,v,vdim);
   result |= CheckIsSubspace(noc, ffGetPtr(wrk2,nor1,noc),nor2,w,wdim);

   sysFree(wrk1);
   sysFree(wrk2);
   sysFree(v);
   sysFree(w);
   sysFree(piv);
   return result;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Matrix_SumIntersection(int q)
{
    int result = 0;
    for (int noc = 1; result == 0 && noc < 100; noc += noc / 10 + 1) {
	result |= TestSumInter1(noc);
	result |= TestSumInter2(noc);
    }
    return result;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
