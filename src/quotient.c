////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Projection on quotient
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup spinup
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Projection on quotient.
/// This function calculates the projection of a matrix onto the quotient by a
/// subspace. The first matrix, @a subspace must be in echelon form, while the
/// second argument can be any matrix. Of course both matrices must be over the
/// same field and have the same number of columns. The return value is a
/// pointer to a matrix containing the projections the @a vectors. This matrix
/// is not in echelon form and may even contain null rows.
///
/// The projection depends on the basis for the subspace and is calculated as
/// follows. Let V=F<sup>n×n</sup> and (w<sub>1</sub>,...w<sub>s</sub>) be a basis
/// for the subspace W≤V. The basis, written as a matrix of row vectors,
/// is assumed to be in semi-echelon form. By looking at the pivot columns we
/// can construct the vectors w<sub>s+1</sub>,...w<sub>n</sub> by taking all vectors which
/// have a exactly one 1 at any non-pivot position and are zero otherwise.
/// Then, (w<sub>1</sub>,...,w<sub>s</sub>,w<sub>s+1</sub>,...,w<sub>n</sub>)
/// is a basis for V in semi-echelon form and defines the decomposition of any
/// vector into subspace and quotient part.
///
/// @param subspace The invariant subspace.
/// @param vectors The vectors to project.
/// @return Projection of @a vectors on the quotient by @a subspace, or 0 on error.

Matrix_t *QProjection(const Matrix_t *subspace, const Matrix_t *vectors)
{
   int i, sdim, qdim;
   int *non_piv;
   Matrix_t *result;
   PTR tmp;

   /* Check the arguments
      ------------------- */
   if (!MatIsValid(subspace) || !MatIsValid(vectors)) {
      return NULL;
   }
   if ((subspace->Field != vectors->Field) || (subspace->Noc != vectors->Noc)) {
      MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
      return NULL;
   }
   if (subspace->PivotTable == NULL) {
      MTX_ERROR1("%E",MTX_ERR_NOTECH);
      return NULL;
   }

   /* Initialize
      ---------- */
   sdim = subspace->Nor;
   qdim = subspace->Noc - sdim;
   result = MatAlloc(subspace->Field,vectors->Nor,qdim);

   /* Calculate the projection
      ------------------------ */
   FfSetNoc(subspace->Noc);
   tmp = FfAlloc(1);
   non_piv = subspace->PivotTable + subspace->Nor;
   for (i = 0; i < vectors->Nor; ++i) {
      int k;
      PTR q = MatGetPtr(result,i);
      FfCopyRow(tmp,MatGetPtr(vectors,i));
      FfCleanRow(tmp,subspace->Data,sdim,subspace->PivotTable);
      for (k = 0; k < qdim; ++k) {
         FfInsert(q,k,FfExtract(tmp,non_piv[k]));
      }
   }
   SysFree(tmp);

   return result;
}


/// Action on Quotient.
/// Given a subspace U≤F<sup>n</sup> and a matrix A∊F<sup>n×n</sup> that maps
/// U into U, this function calculates the action of the matrix on the
/// quotient F<sup>n</sup>/U.
/// As input, the function expects a basis of the subspace in @em subspace,
/// which must be in echelon form, and the matrix operating on the
/// subspace in @a gen. The result is a square matrix with n-dim(U) rows
/// describing the action of A on the quotient in a randomly chosen
/// basis.
///
/// Before calculating the action, QAction() checks if the arguments
/// are valid matrices, and if they are compatible. Both @a subspace
/// and @a gen must be over the same field, have the same number of
/// columns, and @a gen must be square.
///
/// @param subspace Pointer to an invariant subspace.
/// @param gen Matrix operating on the subspace.
/// @return Action of the generator on the quotient, or 0 on error.

Matrix_t *QAction(const Matrix_t *subspace, const Matrix_t *gen)
{
   int k;
   int dim, sdim, qdim;
   int *piv, *non_piv;
   PTR tmp;
   Matrix_t *action;

   /* Check arguments.
      ---------------- */
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

   /* Initialize
      ---------- */
   dim = subspace->Noc;
   sdim = subspace->Nor;
   qdim = dim - sdim;
   if ((action = MatAlloc(subspace->Field,qdim,qdim)) == NULL) {
      return NULL;
   }

   /* Calculate the action on the quotient
      ------------------------------------ */
   FfSetNoc(dim);
   tmp = FfAlloc(1);
   piv = subspace->PivotTable;
   non_piv = piv + subspace->Nor;
   for (k = 0; k < qdim; ++k) {
      int l;
      PTR qx = MatGetPtr(action,k);
      FfCopyRow(tmp,MatGetPtr(gen,non_piv[k]));
      FfCleanRow(tmp,subspace->Data,sdim,piv);
      for (l = 0; l < qdim; ++l) {
         FfInsert(qx,l,FfExtract(tmp,non_piv[l]));
      }
   }
   SysFree(tmp);

   return action;
}


/// @}