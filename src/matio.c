////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Read a matrix from a file
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>

/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads matrix contents from a file end returns the matrix.
/// Note: this function can only be called after a matrix header has been read. 
/// To simply read a matrix from the file, use @ref matRead instead.

Matrix_t* matReadData(MtxFile_t* f)
{
   const uint32_t objectType = mfObjectType(f);
   if (objectType != MTX_TYPE_MATRIX) {
      mtxAbort(MTX_HERE, "%s: bad type 0x%lx, expected 0x%lx (MATRIX)",
         f->name, (unsigned long) objectType, (unsigned long) MTX_TYPE_MATRIX);
   }
   Matrix_t* m = matAlloc(f->header[0], f->header[1], f->header[2]);
   ffReadRows(f, m->data, m->nor, m->noc);

   // Make sure a second read attempt will fail.
   f->header[0] = 0xFFFFFFFF;

   return m;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads a matrix from a file and returns the matrix.
/// The given file must have been be opened for reading, see @ref mfOpen.

Matrix_t* matRead(MtxFile_t* f)
{
   mfReadHeader(f);
   return matReadData(f);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Read a matrix from a named file and returns the matrix.

Matrix_t* matLoad(const char* fn)
{
   MtxFile_t* f = mfOpen(fn, "rb");
   Matrix_t* m = matRead(f);
   mfClose(f);
   return m;
}

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
