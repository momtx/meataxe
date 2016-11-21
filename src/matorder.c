////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Order of a matrix
//
// (C) Copyright 1998-2016 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Order of a matrix.
/// This function calculates the order of a matrix. @em mat must be a
/// non-singular, square matrix.
/// Even if @em mat is non-singular, the function may fail. This happens if
/// the order is greater than 1000000, or if the order on any cyclic
/// subspace is greater than 1000.
/// @param mat Pointer to the matrix.
/// @return The order of @em mat, or -1 on error.

int MatOrder(const Matrix_t *mat)
{
   // Check arguments
   if (!MatIsValid(mat)) {
      return -1;
   }
   if (mat->Nor != mat->Noc) {
      MTX_ERROR1("%E",MTX_ERR_NOTSQUARE);
      return -1;
   }

   FfSetField(mat->Field);
   FfSetNoc(mat->Noc);
   const int nor = mat->Nor;
   PTR const m1 = FfAlloc(nor);
   PTR const basis = FfAlloc(nor + 1);
   int * const piv = NALLOC(int,nor + 1);
   char * const done = NALLOC(char,nor);
   PTR const v1 = FfAlloc(1);
   PTR const v2 = FfAlloc(1);
   PTR const v3 = FfAlloc(1);

   if (m1 == NULL || basis == NULL || piv == NULL || done == NULL
       || v1 == NULL || v2 == NULL || v3 == NULL) {
      if (m1 != NULL) SysFree(m1);
      if (basis != NULL) SysFree(basis);
      if (piv != NULL) SysFree(piv);
      if (done != NULL) SysFree(done);
      if (v1 != NULL) SysFree(v1);
      if (v2 != NULL) SysFree(v2);
      if (v3 != NULL) SysFree(v3);
      return -1;
   }
   PTR bend = basis;
   FEL f1;
   int tord;
   int ord;

   memcpy(m1,mat->Data,FfCurrentRowSize * nor);
   memset(done,0,(size_t)nor);
   tord = ord = 1;
   int dim = 0;
   while (dim < nor && tord <= 1000 && ord <= 1000000) {
      // get next start vector
      int j1;
      for (j1 = 0; j1 < nor && done[j1]; ++j1) {
      }
      if (j1 >= nor) {
         break;                 // done
      }
      FfMulRow(v1,FF_ZERO);
      FfInsert(v1,j1,FF_ONE);

      // calculate order on cyclic subspace
      tord = 0;
      int flag = 1;
      FfCopyRow(v3,v1);
      do {
         FfCopyRow(v2,v3);
         if (flag) {
            FfCopyRow(bend,v3);
            PTR bptr = basis;
            for (int i = 0; i < dim; ++i) {
               f1 = FfExtract(bend,piv[i]);
               if (f1 != 0) {
                  FfAddMulRow(bend,bptr,FfNeg(FfDiv(f1, FfExtract(bptr,piv[i]))));
               }
               FfStepPtr(&bptr);
            }
            if ((piv[dim] = FfFindPivot(bend,&f1)) >= 0) {
               done[piv[dim]] = 1;
               ++dim;
               FfStepPtr(&bend);
            } else {
               flag = 0;
            }
         }
         FfMapRow(v2,m1,nor,v3);
         ++tord;
      }
      while (tord <= 1000 && FfCmpRows(v3,v1) != 0);

      // Order = lcm(orders on cyclic subspaces)
      ord = lcm(ord,tord);
   }

   // Clean up
   SysFree(done);
   SysFree(v1);
   SysFree(v2);
   SysFree(v3);
   SysFree(m1);
   SysFree(basis);
   SysFree(piv);

   if ((tord > 1000) || (ord > 1000000)) {
      return -1;
   }
   return ord;
}


/// @}
