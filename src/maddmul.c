/* ============================= C MeatAxe ==================================
   File:        $Id: maddmul.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Add multiple of a matrix.
   --------------------------------------------------------------------------
   (C) Copyright 1997 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"


/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO 


/// @addtogroup mat
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Add a multiple of a matrix.
/// This function adds a multiple of a matrix to another matrix. The 
/// matrices must be compatible for addition.
/// |MatAddMul()| handles special cases (|coeff| equals 0 or 1) in an
/// intelligent way, so there is no need for the caller to do this.
/// @param dest Matrix to add to.
/// @param src Matrix to add.
/// @param coeff Coefficient.
/// @return The function returns |dest|, or |NULL| on error.

Matrix_t *MatAddMul(Matrix_t *dest, const Matrix_t *src, FEL coeff)
{
    /* Check the aguments
       ------------------ */
    if (!MatIsValid(src) || !MatIsValid(dest))
	return NULL;
    if (dest->Field != src->Field || dest->Nor != src->Nor || 
	dest->Noc != src->Noc)
	return MTX_ERROR1("%E",MTX_ERR_INCOMPAT), NULL;

    /* Handle special cases
       -------------------- */
    if (coeff == FF_ONE)
	MatAdd(dest,src);
    else if (coeff == FF_ZERO)
	;
    else
    {
	/* Add multiple
	   ------------ */
	PTR dp = dest->Data, sp = src->Data;
	int n;
	FfSetField(src->Field);
	FfSetNoc(src->Noc);
	for (n = src->Nor; n > 0; --n)
	{
	    FfAddMulRow(dp,sp,coeff);
	    FfStepPtr(&dp);
	    FfStepPtr(&sp);
	}
    }
    return dest;
}

/// @}
