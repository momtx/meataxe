////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Kronecker product of two matrices
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


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
/// @return The tensor product of @p m1 and @p m2, or NULL on error.

Matrix_t *matTensor(const Matrix_t *m1, const Matrix_t *m2)
{
   Matrix_t *temat;         /* The result */
   PTR x2;                  /* Row pointer for m2 */
   long i2;                 /* Row counter for m2 */
   FEL *rowbuf;             /* Holds one row of <m2> in unpacked form */

   // Check arguments
   matValidate(MTX_HERE, m1);
   matValidate(MTX_HERE, m2);
   MTX_ASSERT(m1->field == m2->field);

   // Allocate the result matrix and workspace
   temat = matAlloc(m1->field,m1->nor * m2->nor,m1->noc * m2->noc);
   if (temat == NULL) {
      return NULL;
   }
   if ((temat->nor == 0) || (temat->noc == 0)) {
      return temat;
   }
   rowbuf = NALLOC(FEL,m2->noc);
   if (rowbuf == NULL) {
      matFree(temat);
      return NULL;
   }

   // Calculate the Kronecker product
   x2 = m2->data;
   for (i2 = 0; i2 < m2->nor; ++i2) {
      int k;
      PTR x1, x3;
      int i1;

      // Unpack the row
      for (k = 0; k < m2->noc; ++k) {
         rowbuf[k] = ffExtract(x2,k);
      }

      // Initialize everything for the inner loop
      x1 = m1->data;
      x3 = matGetPtr(temat,i2);
      MTX_ASSERT(x3 != NULL);

      // Loop through all rows of <m1>
      for (i1 = 0; i1 < m1->nor; ++i1) {
         int k1, k3;

         // Loop through all columns of <m1>
         for (k1 = k3 = 0; k1 < m1->noc; ++k1) {
            int k2;
            FEL f = ffExtract(x1,k1);
            if (f == FF_ZERO) {
               k3 += m2->noc;
            } else if (f == FF_ONE) {
               register FEL *rp = rowbuf;
               for (k2 = 0; k2 < m2->noc; ++k2, ++k3, ++rp) {
                  if (*rp != FF_ZERO) {
                     ffInsert(x3,k3,*rp);
                  }
               }
            } else {
               register FEL *rp = rowbuf;
               for (k2 = 0; k2 < m2->noc; ++k2, ++k3, ++rp) {
                  if (*rp != FF_ZERO) {
                     ffInsert(x3,k3,ffMul(f,*rp));
                  }
               }
            }
         }

         ffStepPtr(&x1, m1->noc);

         x3 = ffGetPtr(x3, m2->nor, temat->noc);
      }

      ffStepPtr(&x2, m2->noc);
   }

   // Clean up
   sysFree(rowbuf);
   return temat;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
