////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Map a vector under a permutation.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"

#ifdef PARANOID
#endif



/// @addtogroup ff
/// @{
 
/// Multiply a Vector by a Permutation.
/// This function multiplies the vector @a row from the right with the permutation @a perm
/// and stores the result in @a result. Multiplication of vectors by permutations is defined
/// as follows: if the permutation maps i to k, then the i-ith mark of the vector is stored
/// in the k-th position of the result. 
///
/// Note: @a result and @a row must not overlap. Otherwise the result is undefined.

void ffPermRow(PTR row, const long *perm, PTR result)
{
    register FEL f;
    register int i;
    register const long *p = (long *)perm;

#ifdef PARANOID
    if (row == result)
	mtxAbort(MTX_HERE,"row = result: undefined result!");
#endif

    for (i = 0; i < ffNoc; ++i)
    {
#ifdef PARANOID
	if (*p < 0 || *p > ffNoc)
	    mtxAbort(MTX_HERE,"Invalid point %d in permutation, noc=%d",(int)*p,ffNoc);
#endif
	f = ffExtract(row,i);
	ffInsert(result,*p++,f);
    }
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
