/* ============================= C MeatAxe ==================================
   File:        $Id: polcmp.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Compare polynomials.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"
#include <string.h>

   
/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO




/// Compare Two Plynomials.
/// This function compares two polynomials and returns 0 if the polynomials are equal,
/// -1 if a<b, or 1 if a>b.  The ordering of polynomials is defined as follows.
/// - If a and b are over different fields, the polynomials over 
///   the larger field is greater. 
/// - Otherwise, if they have different degrees, the polynomial
///   with the higher degree is greater. 
/// - If both field and degree are equal, the result of the comparison is 0
///   if the polynomials are equal. Otherwise it is unspecified if the return value
///   is +1 or -1.
/// @param a First polynomial.
/// @param b Second polynomial.
/// @return 1 if a>b, -1 if a<b, 0 if a=b, or -2 on error.

int PolCompare(const Poly_t *a, const Poly_t *b)

{
    int i;

    /* Check the arguments
       ------------------- */
    if (!PolIsValid(a) || !PolIsValid(b))
    {
	MTX_ERROR1("%E",MTX_ERR_BADARG);
	return -2;
    }

    /* Check if the polynomials are over the same field.
       ------------------------------------------------- */
    if (a->Field > b->Field) 
	return 1;
    if (a->Field < b->Field) 
	return -1;

    /* Check if the polynomials have equal degrees.
       -------------------------------------------- */
    if (a->Degree > b->Degree) 
	return 1;
    if (a->Degree < b->Degree) 
	return -1;

    /* Compare the coefficients
       ------------------------ */
    for (i = a->Degree; i >= 0; --i)
    {
	if (a->Data[i] > b->Data[i]) 
	    return 1;
	if (a->Data[i] < b->Data[i]) 
	    return -1;
    }

    /* The polynomials are equal!
       -------------------------- */
    return 0;
}

/// @}
