////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Print a matrix
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

/*
   MTX_DEFINE_FILE_INFO
 */

/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Print a matrix on stdout.
/// This function prints a matrix on the standard output in readable form.
/// If @a name is not 0, the name followed by an equal sign is printed before the matrix.
/// @param name Name to print before the matrix, or 0.
/// @param m Pointer to the matrix.

void matPrint(const char *name, const Matrix_t *m)
{
   PTR x;
   long i, k;

   matValidate(MTX_HERE, m);
   ffSetField(m->Field);
   ffSetNoc(m->Noc);
   x = m->Data;
   if (name != NULL) { printf("%s=\n",name); }
   for (i = 0; i < m->Nor; ++i) {
      for (k = 0; k < m->Noc; ++k) {
         printf("%d",ffToInt(ffExtract(x,k)));
      }
      printf("\n");
      ffStepPtr(&x, m->Noc);
   }
}


/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
