////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, printing
//
// (C) Copyright 1998-2014 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

/*
   MTX_DEFINE_FILE_INFO
 */

/// @addtogroup bs
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Print a bit string on stdout.
/// This function writes a bit string in readable format, i.e., as a sequence of 0's
/// and 1's, to stdout.
/// @param name Name to print before the matrix, or 0.
/// @param bs The bit string.

void BsPrint(const char *name, const BitString_t *bs)
{
   int i;
   if (name != NULL) { printf("%s=\n",name); }
   for (i = 0; i < bs->Size; ++i) {
      printf("%d",BsTest(bs,i));
   }
   printf("\n");
}


/// @}
