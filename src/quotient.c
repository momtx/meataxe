////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Projection on quotient
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup g_spinup
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
/// @return Projection of @a vectors on the quotient by @a subspace, or NULL on error.

Matrix_t *QProjection(const Matrix_t *subspace, const Matrix_t *vectors)
{
   int i, sdim, qdim;
   uint32_t *non_piv;
   Matrix_t *result;
   PTR tmp;

   // Check the arguments
   matValidate(MTX_HERE, subspace);
   matValidate(MTX_HERE, vectors);
   if ((subspace->field != vectors->field) || (subspace->noc != vectors->noc)) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
   }
   if (subspace->pivotTable == NULL) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTECH);
   }

   // Initialize
   sdim = subspace->nor;
   qdim = subspace->noc - sdim;
   result = matAlloc(subspace->field,vectors->nor,qdim);
   if (result == NULL) {
       return NULL;
   }

   // Calculate the projection
   tmp = ffAlloc(1, subspace->noc);
   if (tmp == NULL) {
       matFree(result);
       return NULL;
   }
   non_piv = subspace->pivotTable + subspace->nor;
   for (i = 0; i < vectors->nor; ++i) {
      int k;
      PTR q = matGetPtr(result,i);
      ffCopyRow(tmp,matGetPtr(vectors,i), subspace->noc);
      ffCleanRow(tmp,subspace->data,sdim,subspace->noc, subspace->pivotTable);
      for (k = 0; k < qdim; ++k) {
         ffInsert(q,k,ffExtract(tmp,non_piv[k]));
      }
   }
   sysFree(tmp);

   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
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

   /* Check arguments.
      ---------------- */
   matValidate(MTX_HERE, subspace);
   matValidate(MTX_HERE, gen);
   if (subspace->noc != gen->nor) {
      mtxAbort(MTX_HERE,"subspace and gen: %s",MTX_ERR_INCOMPAT);
   }
   if (gen->nor != gen->noc) {
      mtxAbort(MTX_HERE,"gen: %s",MTX_ERR_NOTSQUARE);
   }

   /* Initialize
      ---------- */
   dim = subspace->noc;
   sdim = subspace->nor;
   qdim = dim - sdim;
   Matrix_t *action = matAlloc(subspace->field,qdim,qdim);
   if (action == NULL) {
      return NULL;
   }

   /* Calculate the action on the quotient
      ------------------------------------ */
   PTR tmp = ffAlloc(1, dim);
   if (tmp == NULL) {
      matFree(action);
      return NULL;
   }
   const uint32_t* const piv = subspace->pivotTable;
   const uint32_t* const non_piv = piv + subspace->nor;
   for (k = 0; k < qdim; ++k) {
      int l;
      PTR qx = matGetPtr(action,k);
      ffCopyRow(tmp,matGetPtr(gen,non_piv[k]), dim);
      ffCleanRow(tmp,subspace->data,sdim,dim, piv);
      for (l = 0; l < qdim; ++l) {
         ffInsert(qx,l,ffExtract(tmp,non_piv[l]));
      }
   }
   sysFree(tmp);

   return action;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
