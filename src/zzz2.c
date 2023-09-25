////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Basic row operations.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>
#include <stdlib.h>

/// @addtogroup ff
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Characteristic of the current field.
/// Like ffOrder, this variable may be used anywhere, but it must not be modified directly.

int ffChar = 0;

/// The current field order.
/// May be used in expressions but must never modified directly. To change the field,
/// use ffSetField().

uint32_t ffOrder = MTX_NVAL;

/// Field generator.
/// This variable contains a generator for the multiplicative group of the current field.

FEL ffGen = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Allocate row vectors.
/// This function allocates a block of memory for «nor» row vectors of size «noc» over the current
/// field (see «ffSetField()». The rows are initialized with zeroes as described in «ffMulRow()».
/// The memory must be released with «sysFree()» when it is no longer needed. The return value is
/// never NULL, even if «nor» or «noc» is zero.

PTR ffAlloc(int nor, int noc)
{
   register long i;
   MTX_ASSERT(nor >= 0);
   MTX_ASSERT(noc >= 0);

   const size_t rowSize = ffRowSize(noc);
   const size_t req = rowSize * (size_t) nor;

   PTR p = (PTR) sysMalloc(req);

   // Initialize all rows with zeroes.
   char* q = (char *) p;
   for (i = nor; i > 0; --i) {
      ffMulRow((PTR) q, FF_ZERO, noc);
      q += rowSize;
   }
   return p;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Free memory.
/// This function frees a block of memory that has been allocated with
/// ffAlloc() before. If the argument is 0, %ffFree() does nothing.
/// @param x Pointer to memory.

void ffFree(PTR x)
{
   sysFree(x);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Copy a row.
/// This function copies the contents of one row to another row.
/// @param dest Pointer to the destination.
/// @param src Pointer to the source.
/// @param noc Row size (number of colums).

void ffCopyRow(PTR dest, PTR src, int noc)
{
   memcpy(dest, src, ffRowSize(noc));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Exchange two rows
/// This function exchanges the contents of two rows.
/// @param dest Pointer to the first row
/// @param src Pointer to the second row
/// @param noc Row size (number of colums).

void ffSwapRows(PTR dest, PTR src, int noc)
{
   long *p1 = (long *) src;
   long *p2 = (long *) dest;
   for (int i = ffRowSize(noc) / sizeof(long); i > 0; --i) {
      long tmp = *p1;
      *p1++ = *p2;
      *p2++ = tmp;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Get row pointer.
/// This function returns a pointer to a row of a matrix, given the row index.
/// @a base must be a pointer to the beginning of a row, but this need not be the first
/// row of the matrix. For example, <tt>x = ffGetPtr(x,1,noc)</tt> can be used to advance a
/// row pointer to the next row.
///
/// Note: The function does not check if the resulting pointer is still inside the matrix.
/// @see ffStepPtr()
/// @param base Pointer to the first row of the matrix.
/// @param row Row index. The first row has index 0.

PTR ffGetPtr(PTR base, int row, int noc)
{
   return (PTR) ((char *)base + ffSize(row, noc));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Advance to next row.
/// <tt>ffStepPtr(&x, noc)</tt> is equivalent to <tt>x = ffGetPtr(x,1,noc)</tt>.
/// @param x Pointer to the row pointer.
/// @param noc Row size (number of columns).

void ffStepPtr(PTR *x, int noc)
{
   *x = (PTR)((char*)*x + ffRowSize(noc));
}


/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
