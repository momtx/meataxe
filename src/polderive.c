/* ============================= C MeatAxe ==================================
   File:        $Id: polderive.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Polynomial derivation.
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


/** 
 ** @addtogroup poly
 ** @{
 **/


/**
 ** Derive a Polynomial.
 ** This function derives a polynomial. Note that the derived polynomial is
 ** stored in @em pol, replacing the original polynomial. The following piece of
 ** code shows how to keep the original polynomial intact while calculating
 ** the derivative:
 ** @code
 ** Poly_t *pol, *der;
 ** ...
 ** der = PolDerive(PolDup(pol));
 ** @endcode
 ** @param pol Pointer to the polynomial.
 ** @return @em pol.

 **/

Poly_t *PolDerive(Poly_t *pol)

{
    int i, maxdeg = -1;
    register FEL *buf;
    FEL f = FF_ZERO;

    /* Check the argument
       ------------------ */
    if (!PolIsValid(pol))
    {
	MTX_ERROR1("%E",MTX_ERR_BADARG);
	return NULL;
    }

    buf = pol->Data;
    FfSetField(pol->Field);
    for (i = 0; i < pol->Degree; ++i)
    {
	f = FfAdd(f,FF_ONE);
	buf[i] = FfMul(buf[i+1],f);
	if (buf[i] != FF_ZERO) 
	    maxdeg = i;
    }
    pol->Degree = maxdeg;
    return pol;
}
      

/** 
 ** @}
 **/
