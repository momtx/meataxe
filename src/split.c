////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Split a representation
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>

/// @addtogroup g_spinup
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Action on a subspace.
///
/// Given a matrix A∊F<sup>n×n</sup> and a subspace U≤F<sup>n</sup> with UA≤U,
/// this function calculates the action of the matrix on the subspace.
///
/// @p subspace is a basis of U and must be in echelon form.
///
/// @p gen is the matrix. It must have the same number of columns as @p subspace and must
/// operate on U.
/// 
/// The returned matrix is a square matrix with dim(U) rows containing the image of the basis
/// vectors under @p gen, expressed in the given basis.

Matrix_t *subspaceAction(const Matrix_t *subspace, const Matrix_t *gen)
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
   for (uint32_t i = 0; i < subspace->nor; ++i) {
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

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Projection on quotient.
/// Given a subspace U≤V=F<sup>n</sup>, this function computes the projection of arbitrary vectors
/// v∈V on the quotient space, V/U.
///
/// @p subspace is a basis of the subspace and must be in echelon form.
///
/// @p gen is the matrix. It must have the same number of columns as @p subspace, or the function
/// fails. It must also operate on U, but this is not verified.
///
/// The returned matrix contains the projected vectors with respect to a certain basis of V/U,
/// which is uniquely determined by the subspace basis (see below). It has n-dim(U) columns and
/// the same number of rows as @p vectors.
///
/// The result is computed by cleaning the input vectors with the given basis of U and removing
/// all pivot columns from the cleaned vector. In other words, the basis of V/U is
/// (u<sub>1</sub>+U, ..., u<sub>m</sub>+U) where m=n-dim(U) and u<sub>i</sub> is the row vector
/// which has a 1 at the i-th non-pivot column and 0 otherwise.

Matrix_t *quotientProjection(const Matrix_t *subspace, const Matrix_t *vectors)
{

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
   const uint32_t sdim = subspace->nor;
   const uint32_t qdim = subspace->noc - sdim;
   Matrix_t* result = matAlloc(subspace->field,vectors->nor,qdim);
   PTR tmp = ffAlloc(1, subspace->noc);
   const uint32_t * const non_piv = subspace->pivotTable + subspace->nor;

   // Calculate the projection
   for (uint32_t i = 0; i < vectors->nor; ++i) {
      PTR q = matGetPtr(result,i);
      ffCopyRow(tmp,matGetPtr(vectors,i), subspace->noc);
      ffCleanRow(tmp,subspace->data,sdim,subspace->noc, subspace->pivotTable);
      for (uint32_t k = 0; k < qdim; ++k) {
         ffInsert(q,k,ffExtract(tmp,non_piv[k]));
      }
   }
   sysFree(tmp);

   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Action on quotient.
///
/// Given a matrix A∊F<sup>n×n</sup> and a subspace U≤F<sup>n</sup> with UA≤U,
/// this function calculates the action of the matrix on the quotient space, V/U.
///
/// @p subspace is a basis of U and must be in echelon form.
///
/// @p gen is the matrix. It must have the same number of columns as @p subspace, or the function
/// fails. It must also operate on U, but this is not verified.
/// 
/// The result is a square matrix with n-dim(U) rows describing the action of A on the quotient
/// in a basis wich is uniquely determined by @p subspace (see @ref quotientProjection).

Matrix_t *quotientAction(const Matrix_t *subspace, const Matrix_t *gen)
{
   matValidate(MTX_HERE, subspace);
   matValidate(MTX_HERE, gen);
   if (subspace->noc != gen->nor) {
      mtxAbort(MTX_HERE,"subspace and gen: %s",MTX_ERR_INCOMPAT);
   }
   if (gen->nor != gen->noc) {
      mtxAbort(MTX_HERE,"gen: %s",MTX_ERR_NOTSQUARE);
   }
   if (subspace->pivotTable == NULL) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTECH);
   }

   const uint32_t dim = subspace->noc;
   const uint32_t sdim = subspace->nor;
   const uint32_t qdim = dim - sdim;
   Matrix_t *action = matAlloc(subspace->field,qdim,qdim);

   PTR tmp = ffAlloc(1, dim);
   const uint32_t* const piv = subspace->pivotTable;
   const uint32_t* const non_piv = piv + subspace->nor;
   for (uint32_t k = 0; k < qdim; ++k) {
      PTR qx = matGetPtr(action,k);
      ffCopyRow(tmp,matGetPtr(gen,non_piv[k]), dim);
      ffCleanRow(tmp,subspace->data,sdim,dim, piv);
      for (uint32_t l = 0; l < qdim; ++l) {
         ffInsert(qx,l,ffExtract(tmp,non_piv[l]));
      }
   }
   sysFree(tmp);

   return action;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Splits a representation.
/// Given a matrix representation of an algebra A and an A-invariant subspace U, this function
/// calculates two new matrix representations corresponding to the subspace and quotient,
/// respectively.
///
/// @param subspace
///    A basis for the invariant subspace. This matrix must be in echelon form.
/// @param rep
///    The representation.
/// @param sub
///    Pointer to a variable receiving the representation on the subspace.
///    @p sub may be NULL if the subspace representation is not needed. If @p sub is not NULL,
///    @c *sub must be NULL.
/// @param quot
///    Pointer to a variable receiving the representation on the quotient.
///    @p quot may be NULL if the quotient representation is not needed. If @p quot is not NULL,
///    @c *quot must be NULL.
///
/// The function fails if the provided subspace is not invariant under the given representation.
/// However, this check is carried out only if the subspace is calculated, i.e., if @p sub is
/// not NULL. The function also fails if subspace and representation are compatible.
///
/// See also @ref subspaceAction, @ref quotientAction.
///
/// The following example code calculates the submodule generated by a given seed vector. If
/// it is a proper submodule, the representation is split.
/// @code
/// MatRep_t *rep;
/// Matrix_t *seed;
/// Matrix_t *subspace;
/// ...
/// subspace = spinup(seed,rep);
/// if (subspace->nor > 0 && subspace->nor < subspace->noc)
/// {
///     MatRep_t *sub = NULL, *quot = NULL;
///     Split(subspace,rep, &sub, &quot);
/// }
/// @endcode

void split(const Matrix_t* subspace, const MatRep_t* rep, MatRep_t** sub, MatRep_t** quot)
{
   mrValidate(MTX_HERE, rep);
   matValidate(MTX_HERE, subspace);
   MTX_ASSERT(subspace->pivotTable != NULL);

   // Subspace
   if (sub != NULL) {
      MTX_ASSERT(*sub == NULL);
      *sub = mrAlloc(0, NULL, 0);
      for (int g = 0; g < rep->NGen; ++g) {
         Matrix_t* gen = subspaceAction(subspace, rep->Gen[g]);
         mrAddGenerator(*sub, gen, 0);
      }
   }

   // Quotient
   if (quot != NULL) {
      MTX_ASSERT(*quot == NULL);
      *quot = mrAlloc(0, NULL, 0);
      for (int g = 0; g < rep->NGen; ++g) {
         Matrix_t* gen = quotientAction(subspace, rep->Gen[g]);
         mrAddGenerator(*quot, gen, 0);
      }
   }
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
