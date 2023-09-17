////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Matrix null space
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mat
/// @{

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
   nsp = matAlloc(mat->Field,mat->Nor,mat->Nor);
   nsp->PivotTable = NREALLOC(nsp->PivotTable,uint32_t,mat->Nor);

   // calculate the null-space
   dim = znullsp(mat->Data,mat->Nor,mat->Noc,nsp->PivotTable,nsp->Data,flags);
   if (flags) {
      sysFree(nsp->PivotTable);
      nsp->PivotTable = NULL;
   }

   // trim result buffer
   nsp->Nor = dim;
   nsp->Data = (PTR) sysRealloc(nsp->Data,ffRowSize(nsp->Noc) * dim);

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


/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
