////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Action on a subspace.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>
#include <stdlib.h>

MTX_DEFINE_FILE_INFO

/// @addtogroup spinup
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
   if (!MatIsValid(subspace) || !MatIsValid(gen)) {
      return NULL;
   }
   if (subspace->Noc != gen->Nor) {
      MTX_ERROR1("subspace and gen: %E",MTX_ERR_INCOMPAT);
      return NULL;
   }
   if (gen->Nor != gen->Noc) {
      MTX_ERROR1("gen: %E",MTX_ERR_NOTSQUARE);
      return NULL;
   }

   // set up internal variables
   const int dim = subspace->Noc;
   const int sdim = subspace->Nor;
   FfSetField(subspace->Field);
   Matrix_t *action = MatAlloc(FfOrder,sdim,sdim);
   if (action == NULL) {
      return NULL;
   }
   FfSetNoc(dim);
   PTR tmp = FfAlloc(1);
   if (tmp == NULL) {
      MatFree(action);
      return NULL;
   }

   // calculate the action
   for (int i = 0; i < subspace->Nor; ++i) {
      PTR xi = MatGetPtr(subspace,i);
      MTX_FAIL_IF_NOT(xi != NULL);
      PTR yi = MatGetPtr(action,i);
      MTX_FAIL_IF_NOT(yi != NULL);
      FEL f;

      // calculate the image of the <i>-th row of <subspace>
      FfMapRow(xi,gen->Data,dim,tmp);

      // clean the image with the subspace and store coefficients
      int rc = FfCleanRow2(tmp,subspace->Data,sdim,subspace->PivotTable,yi);
      MTX_FAIL_IF_NOT(rc == 0);
      if (FfFindPivot(tmp,&f) >= 0) {
         MTX_ERROR("Split(): Subspace not invariant");
	 SysFree(tmp);
	 MatFree(action);
	 return NULL;
      }
   }

   // Clean up and return the result.
   SysFree(tmp);
   return action;
}


/// @}
