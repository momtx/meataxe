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

static long znullsp(PTR matrix, long nor, int *piv, PTR nsp, int flags)
{
   PTR x, y, a, b;
   int i;
   long noc = ffNoc;
   long dim;
   FEL f;

   // initialize result with identity
   ffSetNoc(nor);
   x = nsp;
   for (i = 0; i < nor; ++i) {
      piv[i] = -1;
      ffMulRow(x,FF_ZERO);
      ffInsert(x,i,FF_ONE);
      ffStepPtr(&x, nor);
   }

   // gaussian elimination
   x = matrix;
   y = nsp;
   for (i = 0; i < nor; ++i) {
      PTR xx = matrix, yy = nsp;
      long k, p;

      for (k = 0; k < i; ++k) {
         ffSetNoc(noc); // not checked since we know noc is valid
         if (((p = piv[k]) >= 0) && ((f = ffExtract(x,p)) != FF_ZERO)) {
            f = ffNeg(ffDiv(f,ffExtract(xx,p)));
            ffSetNoc(noc);
            ffAddMulRow(x,xx,f);
            ffSetNoc(nor);
            ffAddMulRow(y,yy,f);
         }
         ffSetNoc(noc);
         ffStepPtr(&xx, noc);
         ffSetNoc(nor);
         ffStepPtr(&yy, nor);
      }
      ffSetNoc(noc);
      piv[i] = p = ffFindPivot(x,&f);
      ffStepPtr(&x, noc);
      ffSetNoc(nor);
      ffStepPtr(&y, nor);
   }

   // step 2: reduce the null space to echelon form.
   dim = 0;
   x = y = nsp;
   a = b = matrix;
   for (i = 0; i < nor; ++i) {
      if (piv[i] == -1) {
         ffSetNoc(nor);
         if (y != x) { ffCopyRow(y,x); }
         if (flags) {
            ++dim;
         } else {
            ffCleanRow(y,nsp,dim,nor, piv);
            piv[dim++] = ffFindPivot(y,&f);
         }
         ffStepPtr(&y, nor);
      } else {
         ffSetNoc(noc);
         if (b != a) { ffCopyRow(b,a); }
         ffStepPtr(&b, noc);
      }
      ffSetNoc(nor);
      ffStepPtr(&x, nor);
      ffSetNoc(noc);
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

   // allocate workspace
   nsp = matAlloc(mat->Field,mat->Nor,mat->Nor);
   if (nsp == NULL) {
      return NULL;
   }
   nsp->PivotTable = NREALLOC(nsp->PivotTable,int,mat->Nor);
   if (nsp->PivotTable == NULL)
   {
       matFree(nsp);
       return NULL;
   }

   // calculate the null-space
   ffSetNoc(mat->Noc);
   dim = znullsp(mat->Data,mat->Nor,nsp->PivotTable,nsp->Data,flags);
   if (dim == -1)
   {
      matFree(nsp);
      return NULL;
   }
   if (flags) {
      sysFree(nsp->PivotTable);
      nsp->PivotTable = NULL;
   }

   // resize the result buffer to its actual size
   nsp->Nor = dim;
   nsp->Data = (PTR) sysRealloc(nsp->Data,nsp->RowSize * dim);

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
