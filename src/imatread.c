/* ============================= C MeatAxe ==================================
   File:        $Id: imatread.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Read an integer matrix from a file.
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <stdlib.h>

MTX_DEFINE_FILE_INFO


/**
 ** @addtogroup imat
 ** @{
 **/

/**
 ** Read an integer matrix from a file.
 ** This function reads an integer matrix from a file.
 ** @see ImatLoad()
 ** @param f File to read from.
 ** @return Pointer to the matrix, or 0 on error.
 **/

IntMatrix_t *ImatRead(FILE *f)
{
    IntMatrix_t *m;
    long hdr[3];

    if (SysReadLong(f,hdr,3) != 3) 
    {
	MTX_ERROR("Cannot read header");
	return 0;
    }
    if (hdr[0] != -8)	/* HACK: T_IMAT in 2.3 */
    {
	MTX_ERROR("Not an integer matrix");
	return 0;
    }
    m = ImatAlloc(hdr[1],hdr[2]);
    if (m == 0) 
	return 0;
    if (SysReadLong(f,m->Data,m->Nor*m->Noc) != m->Nor * m->Noc)
    {
	ImatFree(m);
	return 0;
    }
    return m;
}


/**
 ** Read an Integer Matrix From a File.
 ** This function opens a file, reads a single integer matrix, and closes 
 ** the file. To read more than one matrix from a file, use ImatRead().
 ** @param fn File name.
 ** @return Pointer to the matrix, or 0 on error.
 **/

IntMatrix_t *ImatLoad(const char *fn)
{
    FILE *f;
    IntMatrix_t *m;

    if ((f = SysFopen(fn,FM_READ)) == 0)
	return 0;
    m = ImatRead(f);
    fclose(f);
    return m;
}

/**
 ** @}
 **/
