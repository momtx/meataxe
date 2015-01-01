////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, core functions
//
// (C) Copyright 1998-2014 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

/// @defgroup bs Bit Strings
/// @{
/// @details
/// The BitString_t type represents a string of 0's and 1's. You may also think of it as a subset
/// {0,1,2,...,N-1}, where each 1 in the bit string means that the subset contains the
/// corresponding element. Unlike the set_t type bit strings are static objects which cannot
/// change their size. However, operations like the intersection of two sets are much faster with
/// BitString_t than with set_t objects.
///
/// Bit strings can be written to and read from files. The file format consists
/// of a header followed by the data. The header contains three 32-bit integers.
/// The first number is $-3$, the second number is the size of the bit
/// string, and the third number is always zero. The file header is followed by
/// the bit string data which is written as a sequence of 32-bit integers where
/// bit 0 of the first integer contains the first bit of the string.

/// @class BitString_t
/// A bit string.

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

#define BS_MAGIC 0x3ff92541

#define BS(size) (((size) + sizeof(long) - 1) / sizeof(long))

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check if a bit string is valid.
/// This function checks if its argument, is a pointer to a valid
/// bit string.
/// If the bit string is o.k., the function returns 1.
/// Otherwise, an error is signaled and, if the error handler does not
/// terminate the program, the function returns 0.
/// @param bs Bit string to check.
/// @return 1 if \p bs is a valid bit string, 0 otherwise.

int BsIsValid(const BitString_t *bs)
{
   if (bs == NULL) {
      MTX_ERROR("NULL bit string");
      return 0;
   }
   if ((bs->Magic != BS_MAGIC) || (bs->Size < 0)) {
      MTX_ERROR2("Invalid bit string (magic=%d, size=%d)",
                 (int)bs->Magic,bs->Size);
      return 0;
   }
   if (bs->BufSize != (int) BS(bs->Size)) {
      MTX_ERROR2("Inconsistent bit string size %d, %d)",
                 bs->Size,(int) BS(bs->Size));
      return 0;
   }
   return 1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Create a bit string.
/// This function creates a new bit string of the specified size. The new bit string is
/// filled with zeros.
/// @see BsFree()
/// @param size Size of the bit string.
/// @return Pointer to the new bit string, or  NULL on error.

BitString_t *BsAlloc(int size)
{
   BitString_t *n;
   int bufsize;         /* Number of long integers in data buffer */

   if (size < 0) {
      MTX_ERROR1("Illegal size %d",size);
      return NULL;
   }

   // allocate data buffer
   bufsize = BS(size);
   n = (BitString_t *) SysMalloc(sizeof(BitString_t) +
                                 (bufsize == 0 ? 0 : bufsize - 1) * sizeof(long));
   if (n == NULL) {
      MTX_ERROR("Cannot allocate bit string");
      return NULL;
   }

   // initialize
   n->Magic = BS_MAGIC;
   n->Size = size;
   n->BufSize = bufsize;
   memset(n->Data,0,sizeof(long) * bufsize);

   return n;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Free a bit string.
/// This function frees a bit string, releasing all associated memory.
/// @see BsAlloc()

int BsFree(BitString_t *bs)
{
   if (!BsIsValid(bs)) {
      return -1;
   }
   memset(bs,0,sizeof(BitString_t));
   SysFree(bs);
   return 0;
}


/// @}