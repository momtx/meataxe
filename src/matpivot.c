////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Rebuild the pivot table of a matrix
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

static int zmkpivot(PTR matrix, int nor, int noc, int *piv, int *ispiv)
{
   PTR x;
   int i, k;

   memset(ispiv,0,sizeof(int) * noc);

   /* Build the pivot table in <piv>.
      Keep track of assigned pivot columns in <ispiv>.
      ------------------------------------------------ */
   for (i = 0, x = matrix; i < nor && i < noc; ++i, FfStepPtr(&x)) {
      FEL f;
      int newpiv = FfFindPivot(x,&f);
      if (ispiv[newpiv]) {
         MTX_ERROR1("%E",MTX_ERR_NOTECH);
         return -1;
      }
      piv[i] = newpiv;
      ispiv[newpiv] = 1;
   }

   /* Insert the non-pivot columns
      ---------------------------- */
   for (k = 0; k < noc; ++k) {
      if (!ispiv[k]) {
         piv[i++] = k;
      }
   }

   MTX_VERIFY(i == noc);
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Reduce to echelon form.
/// This function builds the pivot table of a matrix. Unlike MatEchelonize()
/// this function assumes that @a mat is already in echelon form.
/// @param mat Pointer to the matrix.
/// @return 0 on success, -1 on error.

int MatPivotize(Matrix_t *mat)
{
   int rc;
   int *newtab;
   static int *is_pivot = NULL;
   static int maxnoc = -1;

   // check argument
   if (!MatIsValid(mat)) {
      return -1;
   }

   // (re-)allocate pivot table
   newtab = NREALLOC(mat->PivotTable,int,mat->Noc);
   if (newtab == NULL) {
      MTX_ERROR1("Cannot allocate pivot table (size %d)",mat->Noc);
      return -1;
   }
   mat->PivotTable = newtab;
   if (mat->Noc > maxnoc) {
      int *new_is_piv = NREALLOC(is_pivot,int,mat->Noc);
      if (new_is_piv == NULL) {
         MTX_ERROR("Cannot allocate temporary table");
         return -1;
      }
      is_pivot = new_is_piv;
      maxnoc = mat->Noc;
   }

   // build the pivot table
   FfSetField(mat->Field);
   FfSetNoc(mat->Noc);
   rc = zmkpivot(mat->Data,mat->Nor,mat->Noc,mat->PivotTable,is_pivot);

   return rc;
}


/// @}