////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Matrix nullspace, Gaussian eliminiation, and related functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <string.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static uint32_t zmkechelon(PTR matrix, uint32_t nor, uint32_t noc, uint32_t *piv, uint8_t *ispiv)
{
   uint32_t i;

   // Initialize the table
   for (i = 0; i < noc; ++i) {
      piv[i] = i;
      ispiv[i] = 0;
   }

   // Echelonize the matrix and build the pivot table in <piv>.
   // Keep track of assigned pivot columns in <ispiv>.
   PTR x = matrix;
   uint32_t rank = 0;
   PTR newrow = matrix;
   newrow = matrix;
   for (i = 0, x = matrix; i < nor && rank < noc; ++i, ffStepPtr(&x, noc)) {
      if (rank < i) {
         ffCopyRow(newrow,x, noc);
      }
      ffCleanRow(newrow,matrix,rank,noc,piv);
      FEL f;
      uint32_t newpiv = ffFindPivot(newrow,&f, noc);
      if (newpiv != MTX_NVAL) {
         piv[rank] = newpiv;
         ispiv[newpiv] = 1;
         ++rank;
         ffStepPtr(&newrow, noc);
      }
   }

   // Insert the non-pivot columns
   uint32_t j = rank;
   for (i = 0; i < noc; ++i) {
      if (!ispiv[i]) {
         piv[j++] = i;
      }
   }
   MTX_ASSERT(j == noc);

   return rank;
}


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reduce to echelon form
///
/// This function performs a Gaussian elimination on the matrix @a mat. After return, @a mat is in
/// semi echelon form and its pivot table has been set up.
/// If the rank of @a mat was smaller than the number of rows, some rows are removed during the
/// process. This function can also be used to rebuild the pivot table after the matrix has been
/// modified.
///
/// @param mat Pointer to the matrix.
/// @return Rank (=number of rows) of the echelongized matrix.

uint32_t matEchelonize(Matrix_t *mat)
{
   matValidate(MTX_HERE, mat);

   // Re-allocate the pivot table. This is not really necessary, since
   // «Noc» should never change without releasing the pivot table, but
   // this would be a really nasty bug....
   mat->pivotTable = NREALLOC(mat->pivotTable,uint32_t,mat->noc);

   // Build the pivot table
   uint8_t *is_pivot = NALLOC(uint8_t, mat->noc);
   ffSetField(mat->field);
   const uint32_t rank = zmkechelon(mat->data,mat->nor,mat->noc,mat->pivotTable, is_pivot);
   sysFree(is_pivot);

   // If the rank is less than the number of rows, remove null rows
   if (rank != mat->nor) {
      mat->nor = rank;
      mat->data = (PTR) sysRealloc(mat->data,ffSize(rank, mat->noc));
   }

   return rank;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Nullity of a matrix.
/// This function calculates the dimension of the null-space of a matrix.
/// Unlike matNullity__() this function does not modify the matrix.
/// @param mat Pointer to the matrix.
/// @return Nullity of the matrix, or -1 on error.

uint32_t matNullity(const Matrix_t *mat)
{
   return matNullity__(matDup(mat));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Nullity of a matrix.
/// This function calculates the dimension of the null-space of a matrix
/// and deletes the matrix.
/// @param mat Pointer to the matrix.
/// @return Nullity of @em mat, or -1 on error.

uint32_t matNullity__(Matrix_t *mat)
{
   matEchelonize(mat);
   uint32_t nul = mat->noc - mat->nor;
   matFree(mat);
   return nul;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int zmkpivot(PTR matrix, uint32_t nor, uint32_t noc, uint32_t *piv, uint8_t *ispiv)
{
   // Extract the pivot columns to «piv».
   memset(ispiv,0,sizeof(uint8_t) * noc);
   PTR x = matrix;
   uint32_t i = 0;
   for (i = 0; i < nor && i < noc; ++i, ffStepPtr(&x, noc)) {
      FEL f;
      uint32_t newpiv = ffFindPivot(x,&f, noc);
      if (newpiv == MTX_NVAL || ispiv[newpiv])
         mtxAbort(MTX_HERE, "%s", MTX_ERR_NOTECH);
      piv[i] = newpiv;
      ispiv[newpiv] = 1;
   }

   // Append the non-pivot columns
   for (uint32_t k = 0; k < noc; ++k) {
      if (!ispiv[k]) {
         piv[i++] = k;
      }
   }

   MTX_ASSERT(i == noc);
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Create pivot table.
/// This function creates or updates the pivot table of a matrix. Unlike @ref matEchelonize this
/// function assumes that @a mat is already in echelon form. If it is not, @c matPivotize() fails
/// and aborts the program.

void matPivotize(Matrix_t *mat)
{
   matValidate(MTX_HERE, mat);

   mat->pivotTable = NREALLOC(mat->pivotTable,uint32_t,mat->noc);
   uint8_t *isPivot = NALLOC(uint8_t, mat->noc);
   ffSetField(mat->field);
   zmkpivot(mat->data,mat->nor,mat->noc,mat->pivotTable,isPivot);
   sysFree(isPivot);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Null space.
/// This function calculates the null-space of a matrix. The matrix is passed
/// as first argument, |nor| is the number of rows, |piv| must be a pointer
/// to an array of at least $|nor|+1$ integers, and |nsp| must be a pointer to
/// a square matrix of size |nor|.
///
/// If the function is successfull (non-negative return value),
/// - |matrix| is reduced to echelon form,
/// - |nsp| contains the null-space in echelon form, and
/// - |piv| contains a pivot table for the null space.
/// If |flags| is nonzero, the null-space is not reduced to echelon form,
/// and the contents of |piv| are undefined.
///
/// @return The dimension of the null-space, or -1 on error.

static long znullsp(PTR matrix, int nor, int noc, uint32_t *piv, PTR nsp, int flags)
{
   PTR x, y, a, b;
   FEL f;

   // initialize result with identity
   x = nsp;
   for (uint32_t i = 0; i < nor; ++i) {
      piv[i] = MTX_NVAL;
      ffMulRow(x,FF_ZERO, nor);
      ffInsert(x,i,FF_ONE);
      ffStepPtr(&x, nor);
   }

   // gaussian elimination
   x = matrix;
   y = nsp;
   for (uint32_t i = 0; i < nor; ++i) {
      PTR xx = matrix, yy = nsp;
      uint32_t k, p;

      for (k = 0; k < i; ++k) {
         if (((p = piv[k]) != MTX_NVAL) && ((f = ffExtract(x,p)) != FF_ZERO)) {
            f = ffNeg(ffDiv(f,ffExtract(xx,p)));
            ffAddMulRow(x,xx,f, noc);
            ffAddMulRow(y,yy,f, nor);
         }
         ffStepPtr(&xx, noc);
         ffStepPtr(&yy, nor);
      }
      piv[i] = p = ffFindPivot(x,&f, noc);
      ffStepPtr(&x, noc);
      ffStepPtr(&y, nor);
   }

   // step 2: reduce the null space to echelon form.
   uint32_t dim = 0;
   x = y = nsp;
   a = b = matrix;
   for (uint32_t i = 0; i < nor; ++i) {
      if (piv[i] == MTX_NVAL) {
         if (y != x) 
            ffCopyRow(y,x, nor);
         if (flags) {
            ++dim;
         } else {
            ffCleanRow(y,nsp, dim, nor, piv);
            piv[dim++] = ffFindPivot(y,&f, nor);
         }
         ffStepPtr(&y, nor);
      } else {
         if (b != a) 
            ffCopyRow(b,a, noc);
         ffStepPtr(&b, noc);
      }
      ffStepPtr(&x, nor);
      ffStepPtr(&a, noc);
   }

   return dim;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Null-space of a matrix
/// This function calculates the null-space of a matrix. Unlike matNullSpace(), this function
/// modifies the orginal matrix, but uses less memory since no temporary workspace is allocated.
/// The result is in echelon form.
///
/// @param mat Pointer to the matrix.
/// @param flags If nonzero, the null-space is not reduced to echelon form.
/// @return Pointer to the null-space, or NULL on error.

Matrix_t *matNullSpace_(Matrix_t *mat, int flags)
{
   long dim;
   Matrix_t *nsp;

   // check arguments
   matValidate(MTX_HERE, mat);

   // allocate workspace (sets the field as side effect)
   nsp = matAlloc(mat->field,mat->nor,mat->nor);
   nsp->pivotTable = NREALLOC(nsp->pivotTable,uint32_t,mat->nor);

   // calculate the null-space
   dim = znullsp(mat->data,mat->nor,mat->noc,nsp->pivotTable,nsp->data,flags);
   if (flags) {
      sysFree(nsp->pivotTable);
      nsp->pivotTable = NULL;
   }

   // trim result buffer
   nsp->nor = dim;
   nsp->data = (PTR) sysRealloc(nsp->data,ffRowSize(nsp->noc) * dim);

   return nsp;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Null-space of a matrix
/// This function calculates the null-space of a matrix. Unlike matNullSpace_() and
/// matNullSpace__(), this function does not change the original matrix, but it allocates
/// a temporary copy of the matrix and thus needs more memory.
/// @param mat Pointer to the matrix.
/// @return Pointer to the null-space, or 0 on error.

Matrix_t *matNullSpace(const Matrix_t *mat)
{
   Matrix_t *tmp, *nsp;

   // Check arguments
   matValidate(MTX_HERE, mat);

   // Non-destructive null-space
   if ((tmp = matDup(mat)) == NULL) {
      mtxAbort(MTX_HERE,"Cannot duplicate matrix");
      return NULL;
   }
   nsp = matNullSpace_(tmp,0);
   matFree(tmp);
   return nsp;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Null-space of a matrix
/// This function calculates the null-space of a matrix and deletes the original matrix.
/// @see matNullSpace_(), matNullSpace()
/// @param mat Pointer to the matrix.
/// @return Pointer to the null-space, or 0 on error.

Matrix_t *matNullSpace__(Matrix_t *mat)
{
   Matrix_t *nsp;
   nsp = matNullSpace_(mat,0);
   matFree(mat);
   return nsp;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Clean a matrix.
///
/// This function "cleans" a matrix with a space, i.e., it adds suitable linear combinations
/// of the rows in @em sub to the rows of @em mat such that all pivot columns in @em mat are
/// zero. Both matrices must be over the same field and have the same number of columns.
/// The second matrix, @em sub, must be in echelon form. The cleaned matrix is reduced to
/// echelon form.
/// @return Rank (= number of rows) of the cleaned matrix.

uint32_t matClean(Matrix_t* mat, const Matrix_t* sub)
{
   matValidate(MTX_HERE, mat);
   matValidate(MTX_HERE, sub);
   if (mat->field != sub->field || mat->noc != sub->noc)
      mtxAbort(MTX_HERE, "%s", MTX_ERR_INCOMPAT);
   if (sub->pivotTable == NULL)
      mtxAbort(MTX_HERE, "Subspace: %s", MTX_ERR_NOTECH);

   // Clean
   for (uint32_t i = 0; i < mat->nor; ++i) {
      PTR m = matGetPtr(mat, i);
      ffCleanRow(m, sub->data, sub->nor, sub->noc, sub->pivotTable);
   }

   // Reduce to echelon form
   return matEchelonize(mat);
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
