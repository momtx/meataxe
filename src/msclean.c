////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Matrix set cleaning functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup matset
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

int matFindPivot(const Matrix_t *mat, int *row, int *col, FEL *f)
{
   int i;
   for (i = 0; i < mat->Nor; ++i) {
      FEL g;
      uint32_t piv = ffFindPivot(matGetPtr(mat,i),&g,mat->Noc);
      if (piv != MTX_NVAL) {
         if (f != NULL) { *f = g; }
         *row = i;
         *col = piv;
         return 0;
      }
   }
   return -1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Clean a matrix with a matrix set.
/// This function cleans a matrix with a matrix set by adding suitable
/// multiples of the members of the set to the matrix. When the function
/// returns, all pivot positions of @a mat, as defined by @a set, are zero.
/// @param set Pointer to the matrix set.
/// @param mat Matrix to be cleaned.
/// @return 0 on success, -1 on error.

int msClean(const MatrixSet_t *set, Matrix_t *mat)
{
   int i;
   MatrixSetElement_t *l;

   /* Check the arguments.
      -------------------- */
   matValidate(MTX_HERE, mat);
   if (set->Len > 0) {
      Matrix_t *mat0 = set->List[0].Matrix;
      if ((mat->Field != mat0->Field) || (mat->Nor != mat0->Nor)
          || (mat->Noc != mat0->Noc)) {
         mtxAbort(MTX_HERE,"Cannot clean: %s",MTX_ERR_INCOMPAT);
         return -1;
      }
   }

   /* Clean the matrix.
      ----------------- */
   for (i = 0, l = set->List; i < set->Len; ++i, ++l) {
      FEL f = ffExtract(matGetPtr(mat,l->PivRow),l->PivCol);
      if (f != FF_ZERO) {
         matAddMul(mat,l->Matrix,ffNeg(ffDiv(f,l->PivMark)));
      }
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Extend a matrix set.
/// This function cleans a matrix with a matrix by calling msClean().
/// If the resulting matrix is nonzero, it is added to the matrix set,
/// and the function returns 0. Otherwise the return value is 1,
/// indicating that the matrix is now zero.
///
/// @attention Once a matrix has been added to a matrix set, i.e., after
/// msCleanAndAppend() returns zero, the application is no longer allowed to
/// modify or free the matrix. The matrix will be freed, when the matrix set
/// is destroyed.
/// @param set Pointer to the matrix set.
/// @param mat Matrix to be added.
/// @return 0 if the matrix was added, 1 if the matrix was already in
/// the span of @a set, -1 on error.

int msCleanAndAppend(MatrixSet_t *set, Matrix_t *mat)
{
   int piv_row, piv_col;
   FEL piv_mark;
   MatrixSetElement_t *newlist, *newelem;

   if (msClean(set,mat) != 0) {
      return -1;
   }
   if (matFindPivot(mat,&piv_row,&piv_col,&piv_mark) < 0) {
      return 1;
   }

   /* Extend the matrix list.
      ----------------------- */
   newlist = NREALLOC(set->List,MatrixSetElement_t,set->Len + 1);
   if (newlist == NULL) {
      mtxAbort(MTX_HERE,"Cannot extend matrix set");
      return -1;
   }
   set->List = newlist;
   set->Len++;

   /* Insert the matrix into the list.
      -------------------------------- */
   newelem = set->List + set->Len - 1;
   newelem->Matrix = mat;
   newelem->PivRow = piv_row;
   newelem->PivCol = piv_col;
   newelem->PivMark = piv_mark;

   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
