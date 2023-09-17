////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Read a matrix from a file
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>

/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read matrix contents from a file.
/// @param f File to read from.
/// @param header The object header.
/// @return Pointer to the matrix, or 0 on error.

Matrix_t *matReadData(FILE *f, const uint32_t header[3])
{
   if (header[0] >= MTX_TYPE_BEGIN) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTMATRIX);
   }
   Matrix_t* m = matAlloc(header[0],header[1],header[2]);
   ffReadRows(f,m->Data,m->Nor, m->Noc);
   return m;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read a matrix from a file.
/// @param f File to read from.
/// @return Pointer to the matrix, or 0 on error.

Matrix_t *matRead(FILE *f)
{
   uint32_t header[3];
   sysRead32(f,header,3);
   return matReadData(f, header);
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
