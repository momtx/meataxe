////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check MatInsert()
//
// (C) Copyright 1998-2016 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"

#include <stdlib.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestMatInsert1()
{
   Matrix_t *a;

   Matrix_t* const mat = RndMat(FfOrder,10,10);
   Matrix_t* const zero = MatAlloc(FfOrder,10,10);

   /* Check p(x) = 0 */
   Poly_t* const zeroPoly = PolAlloc(FfOrder,-1);
   a = MatInsert(mat,zeroPoly);
   ASSERT_EQ_INT(MatCompare(a,zero),0);
   MatFree(a);
   a = MatInsert_(MatDup(mat),zeroPoly);
   ASSERT_EQ_INT(MatCompare(a,zero),0);
   MatFree(a);
   PolFree(zeroPoly);

   /* Check p(x) = 1 */
   Matrix_t* const identity = MatId(FfOrder,10);
   Poly_t* const onePoly = PolAlloc(FfOrder,0);
   a = MatInsert(mat,onePoly);
   ASSERT_EQ_INT(MatCompare(a,identity),0);
   MatFree(a);
   a = MatInsert_(MatDup(mat),onePoly);
   ASSERT_EQ_INT(MatCompare(a,identity),0);
   MatFree(a);
   PolFree(onePoly);

   /* Check p(x) = x */
   Poly_t *const idPoly = PolAlloc(FfOrder,1);
   a = MatInsert(mat,idPoly);
   ASSERT_EQ_INT(MatCompare(a,mat),0);
   MatFree(a);
   a = MatInsert_(MatDup(mat),idPoly);
   ASSERT_EQ_INT(MatCompare(a,mat),0);
   MatFree(a);
   PolFree(idPoly);

   MatFree(mat);
   MatFree(identity);
   MatFree(zero);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F InsertMatrixIntoPolynomial()
{
   while (NextField() > 0) {
      TestMatInsert1(FfOrder);
   }
}
