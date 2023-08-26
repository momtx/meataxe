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
/// This function compares two row vectors. As with all row operations, the row 
/// length must have been set before with ffSetNoc(). 
/// The return value is negative if the first row is "less" than the second 
/// row, and it is positive if the first row is "greater" than the second 
/// row. However, the ordering defined by ffCmpRows() depends on the internal
/// representation of finite field elements and can differ between different
/// kernels or between different hardware architectures.
/// @param p1 Pointer to the first matrix.
/// @param p2 Pointer to the second matrix.
/// @return The function returns 0 if the two rows are identical. Otherwise the
/// return value is different from 0.

int ffCmpRows(PTR p1, PTR p2)
{	
    return memcmp(p1,p2,ffTrueRowSize(ffNoc));
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
