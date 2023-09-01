////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Basic integer matrix functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


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
   if ((mat->Magic != IMAT_MAGIC) || mat->Nor < 0 || mat->Noc < 0) {
      mtxAbort(sl ? sl : MTX_HERE,"Invalid matrix (nor=%d, noc=%d)", mat->Nor, mat->Noc);
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
   m->Magic = IMAT_MAGIC;
   m->Nor = nor;
   m->Noc = noc;
   m->Data = NALLOC(int32_t,nor * noc);
   if (m->Data == NULL) {
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
   if (mat->Data != NULL) {
      sysFree(mat->Data);
   }
   memset(mat,0,sizeof(IntMatrix_t));
   sysFree(mat);
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
