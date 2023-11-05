////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Order of a matrix
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


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

int matOrder(const Matrix_t *mat)
{
   // Check arguments
   matValidate(MTX_HERE, mat);
   if (mat->nor != mat->noc) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTSQUARE);
      return -1;
   }

   ffSetField(mat->field);
   const int nor = mat->nor;
   PTR const m1 = ffAlloc(nor, nor);
   PTR const basis = ffAlloc(nor + 1, nor);
   uint32_t* const piv = NALLOC(uint32_t ,nor + 1);
   char* const done = NALLOC(char,nor);
   PTR const v1 = ffAlloc(1, nor);
   PTR const v2 = ffAlloc(1, nor);
   PTR const v3 = ffAlloc(1, nor);

   if (m1 == NULL || basis == NULL || piv == NULL || done == NULL
       || v1 == NULL || v2 == NULL || v3 == NULL) {
      if (m1 != NULL) sysFree(m1);
      if (basis != NULL) sysFree(basis);
      if (piv != NULL) sysFree(piv);
      if (done != NULL) sysFree(done);
      if (v1 != NULL) sysFree(v1);
      if (v2 != NULL) sysFree(v2);
      if (v3 != NULL) sysFree(v3);
      return -1;
   }
   PTR bend = basis;
   FEL f1;
   int tord;
   int ord;

   memcpy(m1,mat->data,ffSize(nor,nor));
   memset(done,0,nor);
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
      ffMulRow(v1,FF_ZERO,nor);
      ffInsert(v1,j1,FF_ONE);

      // calculate order on cyclic subspace
      tord = 0;
      int flag = 1;
      ffCopyRow(v3,v1,nor);
      do {
         ffCopyRow(v2,v3, nor);
         if (flag) {
            ffCopyRow(bend,v3, nor);
            PTR bptr = basis;
            for (int i = 0; i < dim; ++i) {
               f1 = ffExtract(bend,piv[i]);
               if (f1 != FF_ZERO) {
                  ffAddMulRow(bend,bptr,ffNeg(ffDiv(f1, ffExtract(bptr,piv[i]))),nor);
               }
               ffStepPtr(&bptr, nor);
            }
            if ((piv[dim] = ffFindPivot(bend,&f1, nor)) != MTX_NVAL) {
               done[piv[dim]] = 1;
               ++dim;
               ffStepPtr(&bend, nor);
            } else {
               flag = 0;
            }
         }
         ffMapRow(v3,v2,m1,nor,nor);
         ++tord;
      }
      while (tord <= 1000 && ffCmpRows(v3,v1, nor) != 0);

      // Order = lcm(orders on cyclic subspaces)
      ord = lcm(ord,tord);
   }

   // Clean up
   sysFree(done);
   sysFree(v1);
   sysFree(v2);
   sysFree(v3);
   sysFree(m1);
   sysFree(basis);
   sysFree(piv);

   if ((tord > 1000) || (ord > 1000000)) {
      return -1;
   }
   return ord;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
