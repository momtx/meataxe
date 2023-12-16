////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check functions for quotient projection.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult QuotientProjection1(int q)
{
   Matrix_t *sub = MkMat(3,5, 1,0,0,0,0, 0,1,1,0,1, 0,0,0,1,1);
   matEchelonize(sub);

   Matrix_t *m1 = MkMat(5,5, 0,0,0,0,1, 0,0,0,1,0, 0,0,1,0,0, 0,1,0,0,0, 1,0,0,0,0);
   Matrix_t *expectedP1 = MkMat(2,2, 0,1, 1,0);
   Matrix_t* p1 = QProjection(sub,m1);
   matEchelonize(p1);
   ASSERT(matCompare(p1,expectedP1) == 0);
   matFree(p1);
   matFree(expectedP1);
   matFree(m1);

   Matrix_t *m2 = MkMat(5,5, 1,1,1,1,1, 1,0,1,0,1, 1,0,1,1,0, 0,1,0,1,1, 0,1,1,1,0);
   Matrix_t *expectedP2 = MkMat(2,2, 0,-1, 1,0);
   Matrix_t* p2 = QProjection(sub,m2);
   matEchelonize(p2);
   ASSERT(matCompare(p2,expectedP2) == 0);
   matFree(p2);
   matFree(expectedP2);
   matFree(m2);

   matFree(sub);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult QuotientProjection2(int q)
{
   int size;
   for (size = 1; size < 100; size += size / 5 + 1) {
      int k;
      Matrix_t *sub = RndMat(ffOrder,size / 2,size);
      Matrix_t *id = matId(ffOrder,size);
      Matrix_t *quot;
      matEchelonize(sub);
      quot = QProjection(sub,id);
      matEchelonize(quot);
      matFree(id);
      for (k = 0; k < 3; ++k) {
         Matrix_t *vec = RndMat(ffOrder,size * 5,size);
         Matrix_t *proj;
         proj = QProjection(sub,vec);
         matEchelonize(proj);
	 ASSERT(IsSubspace(proj,quot,0));
	 ASSERT(proj->nor <= quot->nor);
         matFree(vec);
         matFree(proj);
      }
      matFree(quot);
      matFree(sub);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult QuotientAction(int q)
{
   Matrix_t *sub = MkMat(3,5, 1,0,0,0,0, 0,1,1,0,1, 0,0,0,1,1);
   Matrix_t *expectedOp1 = MkMat(2,2, 1,0, 0,0);
   Matrix_t *expectedOp2 = MkMat(2,2, 1,-1, 0,-2);

   matEchelonize(sub);

   Matrix_t *m1 = MkMat(5,5, 0,0,0,0,1, 0,0,0,1,0, 0,0,1,0,0, 0,1,0,0,0, 1,0,0,0,0);
   Matrix_t* op1 = QAction(sub,m1);
   ASSERT(op1 != NULL);
   ASSERT(matCompare(op1,expectedOp1) == 0);
   matFree(op1);
   matFree(expectedOp1);
   matFree(m1);

   Matrix_t *m2 = MkMat(5,5, 1,1,1,1,1, 1,0,1,0,1, 1,0,1,1,0, 0,1,0,1,1, 0,1,1,1,0);
   Matrix_t* op2 = QAction(sub,m2);
   ASSERT(op2 != NULL);
   ASSERT(matCompare(op2,expectedOp2) == 0);
   matFree(op2);
   matFree(expectedOp2);
   matFree(m2);

   matFree(sub);
   return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
