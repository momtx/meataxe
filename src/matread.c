/* ============================= C MeatAxe ==================================
   File:        $Id: matread.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Read a matrix from a file.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <stdlib.h>

MTX_DEFINE_FILE_INFO

/// @addtogroup mat
/// @{


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read a matrix from a file.
/// @param f File to read from.
/// @return Pointer to the matrix, or 0 on error.

Matrix_t *MatRead(FILE *f)
{
    Matrix_t *m;
    long hdr[3];

    if (SysReadLong(f,hdr,3) != 3) 
    {
	MTX_ERROR("Cannot read header");
	return NULL;
    }
    if (hdr[0] < 2) 
    {
	MTX_ERROR1("%E",MTX_ERR_NOTMATRIX);
	return NULL;
    }
    m = MatAlloc(hdr[0],hdr[1],hdr[2]);
    if (m == NULL) 
	return NULL;
    if (FfReadRows(f,m->Data,m->Nor) != m->Nor)
    {
	MatFree(m);
	return NULL;
    }
    return m;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read a matrix from a file.
/// This function opens a file, reads a single matrix, and closes the file.
/// To read more than one matrix from a file, use |MatRead()|.
/// @param fn File name.
/// @return Pointer to the matrix, or |NULL| on error.

Matrix_t *MatLoad(const char *fn)
{
    FILE *f;
    Matrix_t *m;

    if ((f = SysFopen(fn,FM_READ)) == NULL)
	return NULL;
    m = MatRead(f);
    fclose(f);
    return m;
}


/// @}
