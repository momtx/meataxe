////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Kronecker product of two matrices
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

/// @defgroup tp Tensor Products
/// @{
/// @details
/// These function are used for calculations with tensor products.

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Tensor Product.
/// This function calculates the (Kronecker) tensor product m1⊗m2.
/// Both matrices must be over the same field.
/// @param m1 Pointer to the first matrix.
/// @param m2 Pointer to the second matrix.
/// @return The tensor product of @a m1 and @a m2, or 0 on error.

Matrix_t *MatTensor(const Matrix_t *m1, const Matrix_t *m2)
{
   Matrix_t *temat;         /* The result */
   PTR x2;                  /* Row pointer for m2 */
   long i2;                 /* Row counter for m2 */
   FEL *rowbuf;             /* Holds one row of <m2> in unpacked form */

   /* Check arguments
      --------------- */
#ifdef DEBUG
   if (!MatIsValid(m1) || !MatIsValid(m2)) {
      return NULL;
   }
#endif
   if (m1->Field != m2->Field) {
      MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
      return NULL;
   }

   /* Allocate the result matrix and workspace
      ---------------------------------------- */
   temat = MatAlloc(m1->Field,m1->Nor * m2->Nor,m1->Noc * m2->Noc);
   if (temat == NULL) {
      return NULL;
   }
   if ((temat->Nor == 0) || (temat->Noc == 0)) {
      return temat;
   }
   rowbuf = NALLOC(FEL,m2->Noc);
   if (rowbuf == NULL) {
      MatFree(temat);
      return NULL;
   }

   /* Calculate the Kronecker product
      ------------------------------- */
   x2 = m2->Data;
   for (i2 = 0; i2 < m2->Nor; ++i2) {
      int k;
      PTR x1, x3;
      int i1;

      /* Unpack the row
         -------------- */
      for (k = 0; k < m2->Noc; ++k) {
         rowbuf[k] = FfExtract(x2,k);
      }

      /* Initialize everything for the inner loop
         ---------------------------------------- */
      x1 = m1->Data;
      x3 = MatGetPtr(temat,i2);
      FfSetNoc(temat->Noc);

      /* Loop through all rows of <m1>
         ----------------------------- */
      for (i1 = 0; i1 < m1->Nor; ++i1) {
         int k1, k3;

         /* Loop through all columns of <m1>
            --------------------------------- */
         for (k1 = k3 = 0; k1 < m1->Noc; ++k1) {
            int k2;
            FEL f = FfExtract(x1,k1);
            if (f == FF_ZERO) {
               k3 += m2->Noc;
            } else if (f == FF_ONE) {
               register FEL *rp = rowbuf;
               for (k2 = 0; k2 < m2->Noc; ++k2, ++k3, ++rp) {
                  if (*rp != FF_ZERO) {
                     FfInsert(x3,k3,*rp);
                  }
               }
            } else {
               register FEL *rp = rowbuf;
               for (k2 = 0; k2 < m2->Noc; ++k2, ++k3, ++rp) {
                  if (*rp != FF_ZERO) {
                     FfInsert(x3,k3,FfMul(f,*rp));
                  }
               }
            }
         }

         /* Next row of <m1>
            ---------------- */
         /*x1 = FfGetPtr(x1,1,m1->Noc);*/
         x1 = (PTR)((char *)x1 + m1->RowSize);
         /*x3 = FfGetPtr(x3,m2->Nor,temat->Noc);*/
         x3 = (PTR)((char *)x3 + m2->Nor * temat->RowSize);
      }

      /* Next row of <m2>
         ---------------- */
      /*x2 = FfGetPtr(x2,1,m2->Noc);*/
      x2 = (PTR)((char *)x2 + m2->RowSize);
   }

   /* Clean up
      -------- */
   FREE(rowbuf);
   return temat;
}


/// @}
