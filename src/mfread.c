/* ============================= C MeatAxe ==================================
   File:        $Id: mfread.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Read rows from a file.
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"
#include <string.h>

   
/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO

/// @addtogroup mf
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read row vectors from a file.
/// This function reads @a nrows rows from a data file into a buffer. 
/// Unlike FfReadRows(), this function changes the current row size to
/// the appropriate value, which is stored in the MtxFile_t object.
/// @param f Pointer to the file.
/// @param buf Data buffer.
/// @param nrows Number of rows to read.
/// @return Number of rows that were actually read. Any value other than count indicates an error.

int MfReadRows(MtxFile_t *f, PTR buf, int nrows)
{
    int i;
    register char *b = (char *) buf;

    if (!MfIsValid(f))
	return -1;
    if (FfNoc != f->Noc)
	FfSetNoc(f->Noc);

    /* Handle special case <FfNoc> = 0
       ------------------------------- */
    if (FfNoc == 0)
	return nrows;

    /* Read rows one by one
       -------------------- */
    for (i = 0; i < nrows; ++i)
    {
        if (fread(b,FfTrueRowSize(FfNoc),1,f->File) != 1) break;
	b += FfCurrentRowSize;
    }
    if (ferror(f->File)) 
	MTX_ERROR1("%s: Read failed: %S",f->Name);
    return i;
}

/// @}

