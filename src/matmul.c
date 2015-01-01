////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Matrix multiplication
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply matrices
/// This function multiplies @em dest from the right by @em src.
/// The matrices must be compatible for multiplication, i.e. they must be over
/// the same field, and the number of columns of @em dest must be equal to the
/// number of rows of @em src.
/// The result of the multiplication is stored in @em dest, overwriting the
/// original contents.
/// @see MatPower()
/// @param dest Left factor and result.
/// @param src Right factor.
/// @return The function returns @em dest, or 0 on error.

Matrix_t *MatMul(Matrix_t *dest, const Matrix_t *src)
{
   PTR x, tmp, result;
   long i;

   // check arguments
#ifdef DEBUG
   if (!MatIsValid(src) || !MatIsValid(dest)) {
      MTX_ERROR1("%E",MTX_ERR_BADARG);
      return NULL;
   }
   if ((src->Field != dest->Field) || (src->Nor != dest->Noc)) {
      MTX_ERROR7("Can't multiply %dx%d/GF(%d) by %dx%d/GF(%d): %E",
                 dest->Nor,dest->Noc,dest->Field,src->Nor,src->Noc,src->Field,
                 MTX_ERR_INCOMPAT);
      return NULL;
   }
#endif

   // matrix multiplication
   FfSetField(src->Field);
   FfSetNoc(src->Noc);
   result = tmp = FfAlloc(dest->Nor);
   if (result == NULL) {
      MTX_ERROR("MatMul(): Cannot allocate product");
      return NULL;
   }
   x = dest->Data;
   for (i = 0; i < dest->Nor; ++i) {
      FfMapRow(x,src->Data,src->Nor,tmp);
      FfStepPtr(&tmp);
/*	x = FfGetPtr(x,1,dest->Noc);*/
      x = (PTR)((char*) x + dest->RowSize);
   }
   SysFree(dest->Data);
   dest->Data = result;
   dest->Noc = src->Noc;
   dest->RowSize = FfCurrentRowSize;

   Mat_DeletePivotTable(dest);

   return dest;
}


/// @}