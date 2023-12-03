////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Integer matrices
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

#define IMAT_MAGIC 0x396AA2F2

/// @defgroup imat Integer matrices
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @class IntMatrix_t
/// An integer matrix.
/// The IntMatrix_t structure represents a matrix with integer entries.
/// Both @c Nor and @c Noc may be zero. In this case, @c Data ist still a valid pointer,
/// but the memory block it points to has size zero.

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Checks an integer matrix and aborts the program if the marix is not valid.

void imatValidate(const struct MtxSourceLocation* sl, const IntMatrix_t *mat)
{
   if (mat == NULL) {
      mtxAbort(sl ? sl : MTX_HERE,"NULL matrix");
   }
   if ((mat->typeId != IMAT_MAGIC) || mat->nor < 0 || mat->noc < 0) {
      mtxAbort(sl ? sl : MTX_HERE,"Invalid matrix (nor=%d, noc=%d)", mat->nor, mat->noc);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Creates a new integer matrix. See also @ref imatFree.
///
/// @param nor Number of rows.
/// @param noc Number of columns.
/// @return Pointer to the new matrix or 0 on error.

IntMatrix_t *imatAlloc(int nor, int noc)
{
   IntMatrix_t *m;

   MTX_ASSERT(nor >= 0);
   MTX_ASSERT(noc >= 0);

   // allocate
   m = ALLOC(IntMatrix_t);
   if (m == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate IntMatrix_t structure");
      return NULL;
   }

   // initialize
   m->typeId = IMAT_MAGIC;
   m->nor = nor;
   m->noc = noc;
   m->data = NALLOC(int32_t,nor * noc);
   if (m->data == NULL) {
      sysFree(m);
      mtxAbort(MTX_HERE,"Cannot allocate matrix data");
      return NULL;
   }
   return m;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Destroys an integer matrix, releasing the associated memory.

void imatFree(IntMatrix_t *mat)
{
   imatValidate(MTX_HERE, mat);
   if (mat->data != NULL) {
      sysFree(mat->data);
   }
   memset(mat,0,sizeof(IntMatrix_t));
   sysFree(mat);
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
