////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Multiply matrix by scalar
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"



/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply a Matrix by a Constant.
///
/// @param dest Pointer to the matrix.
/// @param coeff Value to multiply with.
/// @return The function returns @p dest or NULL on error.

Matrix_t *matMulScalar(Matrix_t *dest, FEL coeff)
{
    matValidate(MTX_HERE, dest);

    if (coeff == FF_ONE)
    {
	/* Nothing to do */
    }
    else
    {
	PTR dp = dest->data;
	int n;
	ffSetField(dest->field);
	for (n = dest->nor; n > 0; --n)
	{
	    ffMulRow(dp, coeff, dest->noc);
	    ffStepPtr(&dp, dest->noc);
	}
    }
    return dest;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
