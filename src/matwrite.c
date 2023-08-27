////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Write a matrix into a file
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a matrix to a file.
/// @see MatSave
/// @param mat Pointer to the matrix.
/// @param f Pointer to the file.
/// @return 0 on success, -1 on error.

int matWrite(const Matrix_t *mat, FILE *f)
{
   long hdr[3];

   matValidate(MTX_HERE, mat);
   hdr[0] = mat->Field;
   hdr[1] = mat->Nor;
   hdr[2] = mat->Noc;
   if (sysWriteLong32(f,hdr,3) != 3) {
      mtxAbort(MTX_HERE,"Cannot write header");
      return -1;
   }
   ffSetField(mat->Field);
   if (ffWriteRows(f,mat->Data,mat->Nor, mat->Noc) != mat->Nor) {
      mtxAbort(MTX_HERE,"Cannot write rows");
      return -1;
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a matrix to a file.
/// This function opens a file, writes a matrix to the file, and closes the
/// file. If a file with the specified name already exists, the old contents
/// of the file are destroyd.
/// To write more than one matrix to a file, use matWrite().
/// @param mat Pointer to the matrix.
/// @param fn File name.
/// @return 0 on success, -1 on error.

int matSave(const Matrix_t *mat, const char *fn)
{
   matValidate(MTX_HERE, mat);
   FILE *f;
   if ((f = sysFopen(fn,"wb")) == NULL) {
      mtxAbort(MTX_HERE,"Cannot open %s: %S",fn);
      return -1;
   }
   const int i = matWrite(mat,f);
   fclose(f);
   if (i != 0) {
      mtxAbort(MTX_HERE,"Cannot write matrix to %s",fn);
      return -1;
   }
   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
