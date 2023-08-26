////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check matInsert()
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestMatInsert1()
{
   Matrix_t* const mat = RndMat(ffOrder,10,10);
   Matrix_t* const zero = matAlloc(ffOrder,10,10);

   /* Check p(x) = 0 */
   Poly_t* const zeroPoly = polAlloc(ffOrder,-1);
   Matrix_t* a = matInsert(mat,zeroPoly);
   ASSERT_EQ_INT(matCompare(a,zero),0);
   matFree(a);

   a = matInsert_(matDup(mat),zeroPoly);
   ASSERT_EQ_INT(matCompare(a,zero),0);
   matFree(a);
   polFree(zeroPoly);

   /* Check p(x) = 1 */
   Matrix_t* const identity = matId(ffOrder,10);
   Poly_t* const onePoly = polAlloc(ffOrder,0);
   a = matInsert(mat,onePoly);
   ASSERT_EQ_INT(matCompare(a,identity),0);
   matFree(a);
   a = matInsert_(matDup(mat),onePoly);
   ASSERT_EQ_INT(matCompare(a,identity),0);
   matFree(a);
   polFree(onePoly);

   /* Check p(x) = x */
   Poly_t *const idPoly = polAlloc(ffOrder,1);
   a = matInsert(mat,idPoly);
   ASSERT_EQ_INT(matCompare(a,mat),0);
   matFree(a);
   a = matInsert_(matDup(mat),idPoly);
   ASSERT_EQ_INT(matCompare(a,mat),0);
   matFree(a);
   polFree(idPoly);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Matrix_InsertIntoPolynomial(int q)
{
   int result = TestMatInsert1(ffOrder);
   return result;
}
