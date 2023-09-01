////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Write an integer matrix into a file
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup imat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Writes an integer matrix to a file.

/// @see ImatSave
/// @param mat Pointer to the matrix.
/// @param f Pointer to the file.
/// @return 0 on success, -1 on error.

void imatWrite(const IntMatrix_t *mat, FILE *f)
{
   imatValidate(MTX_HERE, mat);
   uint32_t hdr[3] = {MTX_TYPE_INTMATRIX, mat->Nor, mat->Noc};
   sysWrite32(f,hdr,3);
   sysWrite32(f, mat->Data, (size_t) mat->Nor * mat->Noc);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes an integer matrix to a file. If a file with te given name exists, its contents are
/// replace wth the matrix.
/// See also @ref imatWrite.
///
/// @param mat Pointer to the matrix.
/// @param file_name File name.

void imatSave(const IntMatrix_t *mat, const char *file_name)
{
   imatValidate(MTX_HERE, mat);
   FILE *f = sysFopen(file_name,"wb");
   imatWrite(mat,f);
   fclose(f);
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
