////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check tensor product.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"
#include "c-tensor.h"
#include "c-matrix.h"

#include <stdlib.h>
#include <stdio.h>

/* --------------------------------------------------------------------------
   TestMatId() - Test MatId()
   -------------------------------------------------------------------------- */

static void TestMatTensor2(int fl, int dim)
{
   int i;

   for (i = 0; i < 10; ++i) {
      int i1, i2, k1, k2;
      int nor1 = MtxRandomInt(dim);
      int nor2 = MtxRandomInt(dim);
      int noc1 = MtxRandomInt(dim);
      int noc2 = MtxRandomInt(dim);
      Matrix_t *m1 = RndMat(fl,nor1,noc1);
      Matrix_t *m2 = RndMat(fl,nor2,noc2);
      Matrix_t *m3 = MatTensor(m1,m2);
      FfSetNoc(m3->Noc + m1->Noc + m2->Noc);    /* Avoid errors */
      for (i1 = 0; i1 < nor1; ++i1) {
         PTR r1 = MatGetPtr(m1,i1);
         for (i2 = 0; i2 < nor2; ++i2) {
            PTR r2 = MatGetPtr(m2,i2);
            PTR r3 = MatGetPtr(m3,i1 * nor2 + i2);
            for (k1 = 0; k1 < noc1; ++k1) {
               FEL f1 = FfExtract(r1,k1);
               for (k2 = 0; k2 < noc2; ++k2) {
                  FEL f2 = FfExtract(r2,k2);
                  FEL f3 = FfExtract(r3,k1 * noc2 + k2);
                  if (FfMul(f1,f2) != f3) {
                     Error("Bad value at (%d,%d)x(%d,%d)",i1,k1,i2,k2);
                  }
               }
            }
         }
      }
      MatFree(m1);
      MatFree(m2);
      MatFree(m3);
   }
}


static void TestMatTensor1(int fl)
{
   int dim;
   for (dim = 1; dim < 50; dim += dim / 3 + 1) {
      TestMatTensor2(fl,dim);
   }
}


void TestMatTensor(unsigned flags)
{
   MtxRandomInit(0);

#if 0
   {
      Matrix_t *a, *b;
      int i;
      int t0, t1;

      SelectField(2);
      a = RndMat(2,50,50);
      b = RndMat(2,50,50);
      t0 = SysTimeUsed();
      for (i = 0; i < 10; ++i) {
         MatFree(MatTensor(a,b));
         printf(".");
         fflush(stdout);
      }
      t1 = SysTimeUsed() - t0;
      printf("GF(2): %d.%d\n",t1 / 10,t1 % 10);
      MatFree(a);
      MatFree(b);

      SelectField(3);
      a = RndMat(3,50,50);
      b = RndMat(3,50,50);
      t0 = SysTimeUsed();
      for (i = 0; i < 10; ++i) {
         MatFree(MatTensor(a,b));
         printf(".");
         fflush(stdout);
      }
      t1 = SysTimeUsed() - t0;
      printf("GF(3): %d.%d\n",t1 / 10,t1 % 10);
      MatFree(a);
      MatFree(b);

      exit(0);

   }
#endif

   while (NextField() > 0) {
      TestMatTensor1(FfOrder);
   }
   flags = 0;
}