/* ============================= C MeatAxe ==================================
   File:        $Id: matid.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Identity matrix.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"


/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO 

/**
 ** @addtogroup mat
 ** @{
 **/

/**
 ** Identity matrix
 ** This function creates an identity matrix with @em nor nows over GF(@em fl).
 ** @param fl Field order.
 ** @param nor Number of rows.
 ** @return Pointer to the matrix, or 0 on error.
 ** @see MatAlloc
 **/

Matrix_t *MatId(int fl, int nor)

{
    Matrix_t *m;
    PTR x;
    long i;

    /* Check the arguments
       ------------------- */
    if (fl < 2 || nor < 0)
    {
	MTX_ERROR3("Matid(%d,%d): %E",fl,nor,MTX_ERR_BADARG);
	return NULL;
    }

    /* Allocate an empty matrix
       ------------------------ */
    m = MatAlloc(fl,nor,nor);
    if (m == NULL) 
	return NULL;

    /* Set diagonal elements to 1
       -------------------------- */
    for (i = 0, x = m->Data; i < nor; ++i, FfStepPtr(&x))
	FfInsert(x,i,FF_ONE);

    return m;
}

/**
 ** @}
 **/
