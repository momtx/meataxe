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
   if (mat->Nor != mat->Noc) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTSQUARE);
      return -1;
   }

   ffSetField(mat->Field);
   ffSetNoc(mat->Noc);
   const int nor = mat->Nor;
   PTR const m1 = ffAlloc(nor, nor);
   PTR const basis = ffAlloc(nor + 1, nor);
   int * const piv = NALLOC(int,nor + 1);
   char * const done = NALLOC(char,nor);
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

   memcpy(m1,mat->Data,ffSize(nor,mat->Noc));
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
      ffMulRow(v1,FF_ZERO);
      ffInsert(v1,j1,FF_ONE);

      // calculate order on cyclic subspace
      tord = 0;
      int flag = 1;
      ffCopyRow(v3,v1);
      do {
         ffCopyRow(v2,v3);
         if (flag) {
            ffCopyRow(bend,v3);
            PTR bptr = basis;
            for (int i = 0; i < dim; ++i) {
               f1 = ffExtract(bend,piv[i]);
               if (f1 != FF_ZERO) {
                  ffAddMulRow(bend,bptr,ffNeg(ffDiv(f1, ffExtract(bptr,piv[i]))));
               }
               MTX_ASSERT(ffNoc == mat->Noc, -1);
               ffStepPtr(&bptr, mat->Noc);
            }
            if ((piv[dim] = ffFindPivot(bend,&f1)) >= 0) {
               done[piv[dim]] = 1;
               ++dim;
               MTX_ASSERT(ffNoc == mat->Noc, -1);
               ffStepPtr(&bend, mat->Noc);
            } else {
               flag = 0;
            }
         }
         ffMapRow(v2,m1,nor,v3);
         ++tord;
      }
      while (tord <= 1000 && ffCmpRows(v3,v1) != 0);

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
