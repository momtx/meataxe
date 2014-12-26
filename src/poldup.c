/* ============================= C MeatAxe ==================================
   File:        $Id: poldup.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Duplicate a polynomial.
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

/// @addtogroup poly
/// @{

/// Duplicate a Polynomial.
/// This function creates a copy of an existing polynomial.
/// @param p Pointer to the polynomial.
/// @return A copy of @em p or 0 on error.
/// @see PolAlloc

Poly_t *PolDup(const Poly_t *p)
{
    Poly_t *y;

    if (!PolIsValid(p))
	return NULL;
    y = PolAlloc(p->Field,p->Degree);
    if (y == NULL)
    {
	MTX_ERROR("Cannot allocate polynomial");
	return NULL;
    }
    memcpy(y->Data,p->Data,(p->Degree + 1) * sizeof(FEL));
    return y;
}

/// @}
