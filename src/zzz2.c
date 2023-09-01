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

int ffOrder = -1;

/// Field generator.
/// This variable contains a generator for the multiplicative group of the current field.

FEL ffGen = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Allocate memory and initialize
/// This function allocates a block of memory for a vector (if @em nrows is 1)
/// or a matrix over the current field. Memory is initialized to zero as described in @ref ffMulRow.
/// @a nrows may be zero, in which case the function returns a memory block of
/// size zero which must be freed using ffFree().
/// @param nor Number of rows (may be zero).
/// @return Pointer to the memory block.

PTR ffAlloc(int nrows, int noc)
{
   register long i;
   MTX_ASSERT(nrows >= 0);
   MTX_ASSERT(noc >= 0);

   const size_t rowSize = ffRowSize(noc);
   const size_t req = rowSize * (size_t) nrows;

   PTR p = (PTR) sysMalloc(req);

   // Initialize all rows with zeroes.
   char* q = (char *) p;
   for (i = nrows; i > 0; --i) {
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
   if (x != NULL) {
      free(x);
   }
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
