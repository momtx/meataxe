////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Convert vector to matrix.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
//   Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup tp
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Convert vector to matrix.
/// This function converts a vector with m=rc entries into a r by c
/// matrix by filling the matrix from top to bottom and left to right with
/// the entries of the vector. The vector is taken as the n-th row of
/// @a vecs. A new matrix is allocated and returned. @a noc is the number of
/// columns of the result, which must be a divisor of the number of columns
/// of @a vecs.
/// @see MatrixToVector
/// @param vecs List of vectors.
/// @param n Number of the vector to convert.
/// @param noc Desired number of columns.
/// @return 0 on success, -1 on error.

Matrix_t *VectorToMatrix(Matrix_t *vecs, int n, int noc)
{
   int i;
   Matrix_t *result;

   /* Check arguments.
      ---------------- */
   if (!MatIsValid(vecs)) {
      MTX_ERROR1("vecs: %E",MTX_ERR_BADARG);
      return NULL;
   }
   if ((noc > vecs->Noc) || (vecs->Noc % noc != 0)) {
      MTX_ERROR3("noc=%d (vec:%d): %E",noc,vecs->Noc,MTX_ERR_BADARG);
      return NULL;
   }

   /* Convert the vector.
      ------------------- */
   result = MatAlloc(vecs->Field,vecs->Noc / noc,noc);
   if (result == NULL) {
      return NULL;
   }
   for (i = 0; i < result->Nor; ++i) {
      if (MatCopyRegion(result,i,0, vecs,n,i * noc,1,noc) != 0) {
         MTX_ERROR("Copy failed");
      }
   }
   return result;
}


/// @}