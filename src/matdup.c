/* ============================= C MeatAxe ==================================
   File:        $Id: matdup.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Reduce a matrix to semi echelon form.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <string.h>

MTX_DEFINE_FILE_INFO


/// @addtogroup mat
/// @{

/// Duplicate a matrix
/// This function creates a copy of an existing matrix. The caller is 
/// responsible for destroying the copy with MatFree() when it is no 
/// longer needed.
/// @return A copy of the source Matrix, or 0 on error.

Matrix_t *MatDup(const Matrix_t *src)
{
    Matrix_t *m;

#ifdef DEBUG
    if (!MatIsValid(src)) 
	return NULL;
#endif
    m = MatAlloc(src->Field,src->Nor,src->Noc);
    if (m == NULL) 
    {
	MTX_ERROR("Cannot allocate matrix");
	return NULL;
    }
    memcpy(m->Data,src->Data,FfCurrentRowSize * src->Nor);
    return m;
}

/// @}
