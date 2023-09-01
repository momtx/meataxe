////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Read a matrix from a file
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read a matrix from a file.
/// @param f File to read from.
/// @return Pointer to the matrix, or 0 on error.

Matrix_t *matRead(FILE *f)
{
   int32_t hdr[3];
   sysRead32(f,hdr,3);
   if (hdr[0] < 2) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTMATRIX);
   }
   Matrix_t* m = matAlloc(hdr[0],hdr[1],hdr[2]);
   ffReadRows(f,m->Data,m->Nor, m->Noc);
   return m;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read a matrix from a file.
/// This function opens a file, reads a single matrix, and closes the file.
/// To read more than one matrix from a file, use |matRead()|.
/// @param fn File name.
/// @return Pointer to the matrix.

Matrix_t *matLoad(const char *fn)
{
   FILE *f;
   Matrix_t *m;

   if ((f = sysFopen(fn,"rb")) == NULL) {
      return NULL;
   }
   m = matRead(f);
   fclose(f);
   return m;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
