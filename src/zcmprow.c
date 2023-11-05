////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Compare rows.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>
#include <stdlib.h>

/// @addtogroup ff
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Compare two rows.
/// This function compares two row vectors, which must have the same size.
/// The return value is negative if the first row is "less" than the second 
/// row, and it is positive if the first row is "greater" than the second 
/// row. However, the ordering defined by ffCmpRows() depends on the internal
/// representation of finite field elements and can differ between different
/// kernels or between different hardware architectures.
/// @param p1 Pointer to the first matrix.
/// @param p2 Pointer to the second matrix.
/// @param noc Row size (number of columns).
/// @return The function returns 0 if the two rows are identical. Otherwise the
/// return value is different from 0.

int ffCmpRows(PTR p1, PTR p2, int noc)
{	
    return memcmp(p1,p2,ffRowSizeUsed(noc));
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
