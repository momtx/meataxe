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

/// Writes a matrix to a file. See also @ref matSave.

void matWrite(const Matrix_t* mat, MtxFile_t* file)
{
   matValidate(MTX_HERE, mat);
   uint32_t hdr[3] = { mat->field, mat->nor, mat->noc };
   mfWrite32(file, hdr, 3);
   ffSetField(mat->field);
   ffWriteRows(file, mat->data, mat->nor, mat->noc);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes a matrix to a named file.
/// If a file with the specified name already exists, its contents are destroyed.
/// To write more than one matrix to a file, use @ref matWrite.
///
/// @param mat Pointer to the matrix.
/// @param fn File name.
/// @return 0 on success, -1 on error.

void matSave(const Matrix_t* mat, const char* fn)
{
   matValidate(MTX_HERE, mat);
   MtxFile_t* file = mfOpen(fn, "wb");
   matWrite(mat, file);
   mfClose(file);
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
