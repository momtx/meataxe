////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Basic row operations.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>
#include <stdlib.h>

//#define LPR (ffCurrentRowSize / sizeof(long))


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

/// Current row size.
/// Used by all low-level row operations. %ffNoc is updated automatically when the row size
/// is changed with ffSetNoc().

int ffNoc = 0;

/// The number of bytes occupied by a single row in memory.
/// Equal to <tt>ffRowSize(ffNoc)</tt> and always a multiple of sizeof(long).
//size_t ffCurrentRowSize = (size_t) -1;

/// The number of bytes occupied by a row in a data file.
/// Equal to <tt>ffTrueRowSize(ffNoc)</tt>. ffCurrentRowSizeIo can be smaller than
/// ffCurrentRowSize because there is no padding in data files.
size_t ffCurrentRowSizeIo = -1;


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
   MTX_ASSERT(nrows >= 0, NULL);
   MTX_ASSERT(noc >= 0, NULL);

   const size_t rowSize = ffRowSize(noc);
   const size_t req = rowSize * (size_t) nrows;

   PTR p = (PTR) sysMalloc(req);
   ffSetNoc(noc); // TODO remove

   // Initialize all rows with zeroes.
   char* q = (char *) p;
   for (i = nrows; i > 0; --i) {
      ffMulRow((PTR) q,FF_ZERO);
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
/// As with all row operations, the row length must have been set before
/// with |ffSetNoc()|.
/// @param dest Pointer to the destination.
/// @param src Pointer to the source.

void ffCopyRow(PTR dest, PTR src)
{
   long *s = (long *) src, *d = (long *) dest, i;
   const int LPR = ffRowSize(ffNoc) / sizeof(long);
   for (i = LPR; i > 0; --i) {
      *d++ = *s++;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Exchange two rows
/// This function exchanges the contents of two rows. As with all row
/// operations, the row length must have been set before with |ffSetNoc()|.
/// @param dest Pointer to the first row
/// @param src Pointer to the second row

void ffSwapRows(PTR dest, PTR src)
{
   long *p1 = (long *) src, *p2 = (long *) dest, i, l;
   const int LPR = ffRowSize(ffNoc) / sizeof(long);
   for (i = LPR; i > 0; --i) {
      l = *p1;
      *p1++ = *p2;
      *p2++ = l;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

PTR ffGetPtr(PTR base, int row, int noc)
{
   MTX_ASSERT(ffNoc == noc, NULL);
   return (PTR) ((char *)base + ffSize(row, noc));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Advance to next row.
/// <tt>ffStepPtr(&x, noc)</tt> is equivalent to <tt>x = ffGetPtr(x,1,noc)</tt>.
/// @param x Pointer to the row pointer.
/// @param noc Row size (number of columns).

void ffStepPtr(PTR *x, int noc)
{
   MTX_ASSERT(ffNoc == noc,);
   *x = (PTR)((char*)*x + ffRowSize(noc));
}


/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
