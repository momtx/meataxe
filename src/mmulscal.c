////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Multiply matrix by scalar
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <meataxe.h>



/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply a Matrix by a Constant.
/// @param dest Pointer to the matrix.
/// @param coeff Value to multiply with.
/// @return The function returns @a dest.

Matrix_t *MatMulScalar(Matrix_t *dest, FEL coeff)
{
#ifdef DEBUG
    if (!MatIsValid(dest))
	return NULL;
#endif

    if (coeff == FF_ONE)
    {
	/* Nothing to do */
    }
    else
    {
	PTR dp = dest->Data;
	int n;
	FfSetField(dest->Field);
	FfSetNoc(dest->Noc);
	for (n = dest->Nor; n > 0; --n)
	{
	    FfMulRow(dp,coeff);
	    FfStepPtr(&dp);
	}
    }
    return dest;
}

/// @}
