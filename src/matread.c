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
   Matrix_t *m;
   long hdr[3];

   if (sysReadLong32(f,hdr,3) != 3) {
      mtxAbort(MTX_HERE,"Cannot read header");
      return NULL;
   }
   if (hdr[0] < 2) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTMATRIX);
      return NULL;
   }
   m = matAlloc(hdr[0],hdr[1],hdr[2]);
   if (m == NULL) {
      return NULL;
   }
   if (ffReadRows(f,m->Data,m->Nor, m->Noc) != m->Nor) {
      mtxAbort(MTX_HERE,"File format error: could not read %d rows", m->Nor);
      matFree(m);
      return NULL;
   }
   return m;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read a matrix from a file.
/// This function opens a file, reads a single matrix, and closes the file.
/// To read more than one matrix from a file, use |matRead()|.
/// @param fn File name.
/// @return Pointer to the matrix, or |NULL| on error.

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
