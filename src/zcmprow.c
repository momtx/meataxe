////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Compare rows.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <meataxe.h>
#include <string.h>
#include <stdlib.h>



/// @addtogroup ff
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Compare two rows.
/// This function compares two row vectors. As with all row operations, the row 
/// length must have been set before with FfSetNoc(). 
/// The return value is negative if the first row is "less" than the second 
/// row, and it is positive if the first row is "greater" than the second 
/// row. However, the ordering defined by FfCmpRows() depends on the internal
/// representation of finite field elements and can differ between dirrerent
/// kernels or between different hardware architectures.
/// @param p1 Pointer to the first matrix.
/// @param p2 Pointer to the second matrix.
/// @return The function returns 0 if the two rows are identical. Otherwise the
/// return value is different from 0.

int FfCmpRows(PTR p1, PTR p2)
{	
    return memcmp(p1,p2,FfTrueRowSize(FfNoc));
}

/// @}
