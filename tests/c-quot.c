////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check functions for quotient projection.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"

#include <stdlib.h>
#include <stdio.h>

/*MTX_DEFINE_FILE_INFO*/

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestQuotProj1()
{
   Matrix_t *sub =
      MkMat(3,5, 1,0,0,0,0, 0,1,1,0,1, 0,0,0,1,1);
   Matrix_t *m1 =
      MkMat(5,5, 0,0,0,0,1, 0,0,0,1,0, 0,0,1,0,0, 0,1,0,0,0, 1,0,0,0,0);
   Matrix_t *m2 =
      MkMat(5,5, 1,1,1,1,1, 1,0,1,0,1, 1,0,1,1,0, 0,1,0,1,1, 0,1,1,1,0);
   Matrix_t *p1 =
      MkMat(2,2, 0,1, 1,0);
   Matrix_t *p2 =
      MkMat(2,2, 0,-1, 1,0);
   Matrix_t *prj[2];

   MatEchelonize(sub);
   FfSetNoc(1);
   prj[0] = QProjection(sub,m1);
   prj[1] = QProjection(sub,m2);
   MatEchelonize(prj[0]);
   MatEchelonize(prj[1]);
   if (MatCompare(prj[0],p1) != 0) {
      Error("p1 different");
   }
   if (MatCompare(prj[1],p2) != 0) {
      Error("p2 different");
   }

   MatFree(sub);
   MatFree(m1);
   MatFree(m2);
   MatFree(p1);
   MatFree(p2);
   MatFree(prj[0]);
   MatFree(prj[1]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestQuotProj2()
{
   int size;
   MtxRandomInit(112);
   for (size = 1; size < 100; size += size / 5 + 1) {
      int k;
      Matrix_t *sub = RndMat(FfOrder,size / 2,size);
      Matrix_t *id = MatId(FfOrder,size);
      Matrix_t *quot;
      MatEchelonize(sub);
      FfSetNoc(1);
      quot = QProjection(sub,id);
      MatEchelonize(quot);
      MatFree(id);
      for (k = 0; k < 3; ++k) {
         Matrix_t *vec = RndMat(FfOrder,size * 5,size);
         Matrix_t *proj;
         FfSetNoc(1);
         proj = QProjection(sub,vec);
         MatEchelonize(proj);
         if (!IsSubspace(proj,quot,0) || (proj->Nor > quot->Nor)) {
            Error("Wrong quotient");
         }
         MatFree(vec);
         MatFree(proj);
      }
      MatFree(quot);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F QuotientProjection()
{
   while (NextField() > 0) {
      TestQuotProj1();
      TestQuotProj2();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestQuotOp1()
{
   Matrix_t *sub =
      MkMat(3,5, 1,0,0,0,0, 0,1,1,0,1, 0,0,0,1,1);
   Matrix_t *m1 =
      MkMat(5,5, 0,0,0,0,1, 0,0,0,1,0, 0,0,1,0,0, 0,1,0,0,0, 1,0,0,0,0);
   Matrix_t *m2 =
      MkMat(5,5, 1,1,1,1,1, 1,0,1,0,1, 1,0,1,1,0, 0,1,0,1,1, 0,1,1,1,0);
   Matrix_t *op1 =
      MkMat(2,2, 1,0, 0,0);
   Matrix_t *op2 =
      MkMat(2,2, 1,-1, 0,-2);
   Matrix_t *op[2];

   MatEchelonize(sub);
   if ((op[0] = QAction(sub,m1)) == NULL) {
      Error("QOperation() failed");
   }
   if ((op[1] = QAction(sub,m2)) == NULL) {
      Error("QOperation() failed");
   }
   if (MatCompare(op[0],op1) != 0) {
      Error("op1 different");
   }
   if (MatCompare(op[1],op2) != 0) {
      Error("op2 different");
   }

   MatFree(sub);
   MatFree(m1);
   MatFree(m2);
   MatFree(op1);
   MatFree(op2);
   MatFree(op[0]);
   MatFree(op[1]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F QuotientAction()
{
   while (NextField() > 0) {
      TestQuotOp1();
   }
}
