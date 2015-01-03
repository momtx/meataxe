////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Matrix set cleaning functions
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup matset
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

int MatFindPivot(const Matrix_t *mat, int *row, int *col, FEL *f)
{
   int i;
   for (i = 0; i < mat->Nor; ++i) {
      int piv;
      FEL g;
      piv = FfFindPivot(MatGetPtr(mat,i),&g);
      if (piv >= 0) {
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

int MsClean(const MatrixSet_t *set, Matrix_t *mat)
{
   int i;
   MatrixSetElement_t *l;

   /* Check the arguments.
      -------------------- */
   if (!MsIsValid(set) || !MatIsValid(mat)) {
      return -1;
   }
   if (set->Len > 0) {
      Matrix_t *mat0 = set->List[0].Matrix;
      if ((mat->Field != mat0->Field) || (mat->Nor != mat0->Nor)
          || (mat->Noc != mat0->Noc)) {
         MTX_ERROR1("Cannot clean: %E",MTX_ERR_INCOMPAT);
         return -1;
      }
   }

   /* Clean the matrix.
      ----------------- */
   for (i = 0, l = set->List; i < set->Len; ++i, ++l) {
      FEL f = FfExtract(MatGetPtr(mat,l->PivRow),l->PivCol);
      if (f != FF_ZERO) {
         MatAddMul(mat,l->Matrix,FfNeg(FfDiv(f,l->PivMark)));
      }
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Extend a matrix set.
/// This function cleans a matrix with a matrix by calling MsClean().
/// If the resulting matrix is nonzero, it is added to the matrix set,
/// and the function returns 0. Otherwise the return value is 1,
/// indicating that the matrix is now zero.
///
/// @attention Once a matrix has been added to a matrix set, i.e., after
/// MsCleanAndAppend() returns zero, the application is no longer allowed to
/// modify or free the matrix. The matrix will be freed, when the matrix set
/// is destroyed.
/// @param set Pointer to the matrix set.
/// @param mat Matrix to be added.
/// @return 0 if the matrix was added, 1 if the matrix was already in
/// the span of @a set, -1 on error.

int MsCleanAndAppend(MatrixSet_t *set, Matrix_t *mat)
{
   int piv_row, piv_col;
   FEL piv_mark;
   MatrixSetElement_t *newlist, *newelem;

   if (MsClean(set,mat) != 0) {
      return -1;
   }
   if (MatFindPivot(mat,&piv_row,&piv_col,&piv_mark) < 0) {
      return 1;
   }

   /* Extend the matrix list.
      ----------------------- */
   newlist = NREALLOC(set->List,MatrixSetElement_t,set->Len + 1);
   if (newlist == NULL) {
      MTX_ERROR("Cannot extend matrix set");
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