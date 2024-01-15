////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Action on a subspace.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>


/// @addtogroup g_spinup
/// @{

/// Action on a subspace.
/// Given a subspace U≤F<sup>n</sup> and a matrix A∊F<sup>n×n</sup> that maps
/// U into U, this function calculates the action of the matrix on the subspace.
/// As input, the function expects a basis of the subspace in @a subspace,
/// which must be in chelon form, and the matrix operating on the
/// subspace in @a gen. The result is a square matrix with dim(U) rows
/// containing the image of the basis vectors under A, expressed in the
/// given basis.
///
/// Before calculating the action, SAction() checks if the arguments
/// are valid matrices, and if they are compatible. Both @a subspace
/// and @a gen must be over the same field, have the same number of
/// columns, and @a gen must be square.
///
/// @param subspace Pointer to an invariant subspace.
/// @param gen Matrix operating on the subspace.
/// @return Action of the generator on the subspace, or 0 on error.

Matrix_t *SAction(const Matrix_t *subspace, const Matrix_t *gen)
{
   // check arguments
   matValidate(MTX_HERE, subspace);
   matValidate(MTX_HERE, gen);
   if (subspace->noc != gen->nor) {
      mtxAbort(MTX_HERE,"subspace and gen: %s",MTX_ERR_INCOMPAT);
      return NULL;
   }
   if (gen->nor != gen->noc) {
      mtxAbort(MTX_HERE,"gen: %s",MTX_ERR_NOTSQUARE);
      return NULL;
   }

   // set up internal variables
   const int dim = subspace->noc;
   const int sdim = subspace->nor;
   ffSetField(subspace->field);
   Matrix_t *action = matAlloc(ffOrder,sdim,sdim);
   if (action == NULL) {
      return NULL;
   }
   PTR tmp = ffAlloc(1, dim);
   if (tmp == NULL) {
      matFree(action);
      return NULL;
   }

   // calculate the action
   for (int i = 0; i < subspace->nor; ++i) {
      PTR xi = matGetPtr(subspace,i);
      MTX_ASSERT(xi != NULL);
      PTR yi = matGetPtr(action,i);
      MTX_ASSERT(yi != NULL);
      FEL f;

      // calculate the image of the <i>-th row of <subspace>
      ffMapRow(tmp, xi,gen->data,dim,dim);

      // clean the image with the subspace and store coefficients
      int rc = ffCleanRow2(tmp,subspace->data,sdim,dim,subspace->pivotTable,yi);
      MTX_ASSERT(rc == 0);
      if (ffFindPivot(tmp,&f,dim) != MTX_NVAL) {
         mtxAbort(MTX_HERE,"Split(): Subspace not invariant");
	 sysFree(tmp);
	 matFree(action);
	 return NULL;
      }
   }

   // Clean up and return the result.
   sysFree(tmp);
   return action;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
