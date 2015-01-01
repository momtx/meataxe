////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Print a factored polynomial
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

   
/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO


/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply factored polynomials.
/// Multiplies @em dest by @em src. The previous content of @em dest is lost.
/// @see FpMulP()
/// @param dest Factored polynomial to modify.
/// @param src Factored polynomial.
/// @return The function returns |dest| or |NULL| on error.

FPoly_t *FpMul(FPoly_t *dest, const FPoly_t *src)
{
    int i;

    /* Check the arguments
       ------------------- */
    if (!FpIsValid(src) || !FpIsValid(dest))
	return NULL;

    for (i = 0; i < src->NFactors; ++i)
    {
	if (FpMulP(dest,src->Factor[i],src->Mult[i]) == NULL)
	{
	    MTX_ERROR("Cannot multiply");
	    return NULL;
	}
    }
    return dest;
}


/// @}

