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

void matWrite(const Matrix_t *mat, FILE *f)
{
   matValidate(MTX_HERE, mat);
   uint32_t hdr[3] = {mat->field, mat->nor, mat->noc};
   sysWrite32(f,hdr,3);
   ffSetField(mat->field);
   ffWriteRows(f, mat->data, mat->nor, mat->noc);
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

void matSave(const Matrix_t *mat, const char *fn)
{
   matValidate(MTX_HERE, mat);
   FILE *f = sysFopen(fn,"wb");
   matWrite(mat,f);
   fclose(f);
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
