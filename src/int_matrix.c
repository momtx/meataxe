////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Integer matrices
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>


/// @defgroup imat Integer matrices
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// @class IntMatrix_t
/// A matrix over ℤ, using with 32 bit signed integers.
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

/// Creates an independent copy of an integer matrix.

IntMatrix_t *imatDup(const IntMatrix_t* mat)
{
   imatValidate(MTX_HERE, mat);
   IntMatrix_t *m = imatAlloc(mat->nor, mat->noc);
   memcpy(m->data, mat->data, sizeof(int32_t) * mat->nor * mat->noc);
   return m;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Creates an integer matrix from a given row buffer.
/// The passed row buffer may contain more than @p nor rows. In this case the buffer is resized
/// to the given number of rows.
/// After return the buffer is owned by the matrix and must not be modified except by using the
/// using the imatXxx() functions.

IntMatrix_t* imatCreateFromBuffer(int32_t* buffer, uint32_t nor, uint32_t noc)
{
   IntMatrix_t *m = (IntMatrix_t*) mmAlloc(MTX_TYPE_INTMATRIX, sizeof(IntMatrix_t));
   m->nor = nor;
   m->noc = noc;
   m->data = NREALLOC(buffer, int32_t, sizeof(int32_t) * nor * noc);
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

/// Compares two integer matrices.
/// Returns 0 if the matrices are equal or an nonzero value otherwise.
/// The matrices need not have the same number of rows or columns. Matrices of different dimensions
/// are never equal, though.

int imatCompare(const IntMatrix_t* a, const IntMatrix_t* b)
{
   imatValidate(MTX_HERE, a);
   imatValidate(MTX_HERE, b);
   if (a->nor < b->nor) return -1;
   if (a->nor > b->nor) return 1;
   if (a->noc < b->noc) return -1;
   if (a->noc > b->noc) return 1;
   int32_t* pa = a->data;
   int32_t* pb = b->data;
   for (size_t n = (size_t) a->nor * a->noc; n > 0; --n) {
      if (*pa < *pb) return -1;
      if (*pa > *pb) return 1;
      ++pa;
      ++pb;
   }
   return 0;
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
