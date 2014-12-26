/* ============================= C MeatAxe ==================================
   File:        $Id: mattr.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Transpose a matrix.
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <stdlib.h>

MTX_DEFINE_FILE_INFO

/// @addtogroup mat
/// @{

/// Transpose a matrix.
/// @param src Pointer to the matrix.
/// @return Pointer to the transposed matrix or 0 on error.

Matrix_t *MatTransposed(const Matrix_t *src)
{
    PTR s, d;
    long i;
    Matrix_t *dest;

#ifdef DEBUG
    if (!MatIsValid(src))
    {
	MTX_ERROR1("src: %E",MTX_ERR_BADARG);
	return NULL;
    }
#endif
    dest = MatAlloc(src->Field,src->Noc,src->Nor);
    if (dest == NULL) 
    {
	MTX_ERROR("Cannot allocate result");
	return NULL;
    }
    d = dest->Data;
    for (i = 0; i < src->Noc; ++i)
    {
	int k;
	s = src->Data;
	for (k = 0; k < src->Nor; ++k)
	{
#if defined(DEBUG) && defined(PARANOID)
	    FEL f;
	    FfSetNoc(src->Noc);
	    f = FfExtract(s,i);
	    FfSetNoc(src->Nor);
	    FfInsert(d,k,f);
#else
	    FfInsert(d,k,FfExtract(s,i));
#endif
	    s = (PTR)((char*) s + src->RowSize);
	}
	/*d = FfGetPtr(d,1,dest->Noc);*/
	d = (PTR)((char*) d + dest->RowSize);

    }
    return dest;
}

/// @}
