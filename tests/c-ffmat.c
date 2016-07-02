////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for matrix functions
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMapRow1(PTR mat, PTR a, PTR b, int size)
{
   int i;

   /* Make <mat> the unit matrix.
      --------------------------- */
   for (i = 0; i < size; ++i) {
      PTR m = (PTR)((char *)mat + i * FfCurrentRowSize);
      FfInsert(m,i,FF_ONE);
   }

   for (i = 0; i < size; ++i) {
      int k;

      /* Check that the i-th basis vector is mapped on itself.
         ----------------------------------------------------- */
      FfMulRow(a,FF_ZERO);
      FfInsert(a,i,FF_ONE);
      FfMapRow(a,mat,size,b);
      for (k = 0; k < size; ++k) {
         if ((FfExtract(b,k) != FF_ZERO) ^ (k == i)) {
            TST_FAIL("Error 1");
         }
      }

      FfMapRow(a,mat,i,b);
      for (k = 0; k < size; ++k) {
         if (FfExtract(b,k) != FF_ZERO) {
            TST_FAIL("Error 2");
         }
      }
   }

   for (i = 0; i < size; ++i) {
      FfInsert(a,i,FTab[i % FfOrder]);
   }
   FfMapRow(a,mat,size,b);
   for (i = 0; i < size; ++i) {
      ASSERT_EQ_INT(FfExtract(b,i), FTab[i % FfOrder]);
   }

}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F MapRow()
{
   int size = 10;

   while (NextField() > 0) {
      PTR mat, a, b;
      FfSetNoc(size);
      mat = FfAlloc(size);
      a = FfAlloc(1);
      b = FfAlloc(1);
      TestMapRow1(mat,a,b,size);
      SysFree(mat);
      SysFree(a);
      SysFree(b);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void printrows(PTR x, int n)
{
   for (; n > 0; --n) {
      int k;
      for (k = 0; k < FfNoc; ++k) {
         printf("%d",FfExtract(x,k));
      }
      printf("\n");
      FfStepPtr(&x);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestSumInter1()
{
   PTR w1, w2;
   int *piv;
   int nor1, nor2;
   int i, k;

   nor1 = nor2 = FfNoc;
   w1 = FfAlloc(nor1 + nor2);
   w2 = FfAlloc(nor1 + nor2);
   piv = NALLOC(int,nor1 + nor2);
   k = 0;
   for (i = 0; i < nor1; ++i) {
      PTR x = FfGetPtr(w1,i);
      FfInsert(x,k,FF_ONE);
      if (k % 3 == 0) {
         ++k;
      } else {
         k += 2;
      }
      if (k >= FfNoc) {
         k = 0;
      }
   }
   k = 0;
   for (i = nor1; i < nor1 + nor2; ++i) {
      PTR x = FfGetPtr(w1,i);
      FfInsert(x,k,FF_ONE);
      if (k % 3 == 0) {
         k += 2;
      } else {
         ++k;
      }
      if (k >= FfNoc) {
         k = 0;
      }
   }

   FfSumAndIntersection(w1,&nor1,&nor2,w2,piv);
   ASSERT_EQ_INT(nor1, FfNoc);	// whole space
   ASSERT_EQ_INT(nor2, (FfNoc - 1) / 3 + 1);
   for (i = 0; i < nor2; ++i) {
      ASSERT_EQ_INT(piv[nor1 + i] % 3, 0);
   }

   SysFree(w1);
   SysFree(w2);
   SysFree(piv);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void CheckSubspace(PTR u, int udim, PTR v, int vdim)
{
   int *piv = NALLOC(int,vdim);
   int i, vrank;
   PTR x, y, row;

   /* Echelonize v.
      ------------- */
   y = v;
   vrank = 0;
   for (i = 0, x = v; i < vdim; ++i, FfStepPtr(&x)) {
      FEL f;
      int p;
      FfCleanRow(x,v,vrank,piv);
      p = FfFindPivot(x,&f);
      if (p < 0) {
         continue;
      }
      if (i > vrank) {
         FfCopyRow(y,x);
      }
      piv[vrank++] = p;
      FfStepPtr(&y);
   }

   /* Clean u with v.
      --------------- */
   row = FfAlloc(1);
   for (i = 0, x = u; i < udim; ++i, FfStepPtr(&x)) {
      FEL f;
      FfCopyRow(row,x);
      FfCleanRow(row,v,vrank,piv);
      if (FfFindPivot(row,&f) >= 0) {
         printf("Error at %d/%d, Noc=%d\n",i,udim,FfNoc);
         printf("v: (nor=%d, rank=%d)\n",vdim,vrank);
         printrows(v,vdim);
         printf("u: (nor=%d)\n",udim);
         printrows(u,udim);
         TST_FAIL("CheckSubspace() failed");
         break;
      }
   }
   SysFree(row);
   SysFree(piv);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestSumInter2()
{
   PTR v, w, wrk1, wrk2, x;
   int *piv;
   int nor1, nor2;
   int vdim, wdim;
   int i, k;

   /* Allocate buffers.
      ----------------- */
   vdim = nor1 = MtxRandomInt(FfNoc + 1);
   wdim = nor2 = MtxRandomInt(FfNoc + 1);
   v = FfAlloc(vdim);
   w = FfAlloc(wdim);
   wrk1 = FfAlloc(nor1 + nor2);
   wrk2 = FfAlloc(nor1 + nor2);
   piv = NALLOC(int,nor1 + nor2);

   /* Fill with random values.
      ------------------------ */
   for (i = 0, x = v; i < nor1; ++i, FfStepPtr(&x)) {
      for (k = 0; k < FfNoc; ++k) {
         FfInsert(x,k,FTab[MtxRandomInt(FfOrder)]);
      }
   }
   for (i = 0, x = w; i < nor2; ++i, FfStepPtr(&x)) {
      for (k = 0; k < FfNoc; ++k) {
         FfInsert(x,k,FTab[MtxRandomInt(FfOrder)]);
      }
   }
   memcpy(wrk1,v,nor1 * FfCurrentRowSize);
   memcpy(FfGetPtr(wrk1,nor1),w,nor2 * FfCurrentRowSize);

   FfSumAndIntersection(wrk1,&nor1,&nor2,wrk2,piv);

   /*
      printf("v:\n");
      printrows(v,vdim);
      printf("w:\n");
      printrows(w,wdim);
      printf("sum:\n");
      printrows(wrk1,nor1);
      printf("int:\n");
      printrows(FfGetPtr(wrk2,nor1),nor2);
    */

   /* Check relations between V, W, V+w, and VW.
      ------------------------------------------ */
   CheckSubspace(v,vdim,wrk1,nor1);
   CheckSubspace(w,wdim,wrk1,nor1);
   CheckSubspace(FfGetPtr(wrk2,nor1),nor2,v,vdim);
   CheckSubspace(FfGetPtr(wrk2,nor1),nor2,w,wdim);

   SysFree(wrk1);
   SysFree(wrk2);
   SysFree(v);
   SysFree(w);
   SysFree(piv);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

test_F SumIntersection()
{
   while (NextField() > 0) {
      int size;

      for (size = 1; size < 100; size += size / 10 + 1) {
         FfSetNoc(size);
         TestSumInter1();
         TestSumInter2();
      }
   }
}
