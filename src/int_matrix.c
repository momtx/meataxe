////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Integer matrices
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>


/// @defgroup imat Integer matrices
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// @class IntMatrix_t
/// A matrix with (32 bit singed) integer entries.
/// Both @c nor and @c noc may be zero.

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Checks an integer matrix and aborts the program if the marix is not valid.

void imatValidate(const struct MtxSourceLocation* sl, const IntMatrix_t* mat)
{
   if (mat == NULL) {
      mtxAbort(sl ? sl : MTX_HERE, "NULL matrix");
   }
   if ((mat->typeId != MTX_TYPE_INTMATRIX) || mat->nor < 0 || mat->noc < 0) {
      mtxAbort(sl ? sl : MTX_HERE, "Invalid matrix (nor=%lu, noc=%lu)",
         (unsigned long)mat->nor, (unsigned long)mat->noc);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Creates a new integer matrix. See also @ref imatFree.
///
/// @param nor Number of rows.
/// @param noc Number of columns.

IntMatrix_t *imatAlloc(uint32_t nor, uint32_t noc)
{
   IntMatrix_t *m = (IntMatrix_t*) mmAlloc(MTX_TYPE_INTMATRIX, sizeof(IntMatrix_t));
   m->nor = nor;
   m->noc = noc;
   m->data = NALLOC(int32_t,(size_t) nor * noc);
   return m;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Destroys an integer matrix and releases the associated memory.

void imatFree(IntMatrix_t *mat)
{
   imatValidate(MTX_HERE, mat);
   sysFree(mat->data);
   mat->data = NULL;
   mmFree(mat, MTX_TYPE_INTMATRIX);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads an integer matrix from a file. See also @ref imatLoad.

IntMatrix_t* imatRead(MtxFile_t* file)
{
   mfReadHeader(file);
   if (mfObjectType(file) != MTX_TYPE_INTMATRIX) {
      mtxAbort(MTX_HERE, "%s: unexpected object type 0x%lx (expected integer matrix)",
               file->name, (unsigned long) file->header[0]);
   }
   IntMatrix_t* m = imatAlloc(file->header[1], file->header[2]);
   mfRead32(file, m->data, m->nor * m->noc);
   return m;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Read an integer matrix from a file.
/// This function opens a file, reads a single integer matrix, and closes the file.
/// See also use @ref imatRead.
/// @param fn File name.
/// @return Pointer to the matrix.

IntMatrix_t *imatLoad(const char *fn)
{
   MtxFile_t* file = mfOpen(fn, "rb");
   IntMatrix_t* m = imatRead(file);
   mfClose(file);
   return m;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes an integer matrix to a file.
/// See also @ref imatSave.

void imatWrite(const IntMatrix_t *mat, MtxFile_t *f)
{
   imatValidate(MTX_HERE, mat);
   uint32_t hdr[3] = {MTX_TYPE_INTMATRIX, mat->nor, mat->noc};
   mfWrite32(f,hdr,3);
   mfWrite32(f, mat->data, (size_t) mat->nor * mat->noc);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes an integer matrix to a named file. If a file with the same name exists, its contents are
/// replaced with the matrix.
/// See also @ref imatWrite.

void imatSave(const IntMatrix_t* mat, const char* file_name)
{
   imatValidate(MTX_HERE, mat);
   MtxFile_t* f = mfOpen(file_name, "wb");
   imatWrite(mat, f);
   mfClose(f);
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
