/* ============================= C MeatAxe ==================================
   File:        $Id: mfwritelong.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Write long integers.
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




/**
!section obj.file
 ** Write long integers to a file.
 ** @param f
    Pointer to the file.
 ** @param buf
    Data buffer.
 ** @param count
    Number of integers to write.
 ** @return
    Number of integers that were actually written. Any value other than |count|
    indicates an error.
!description
    This function writes |count| long integers from buffer into a data file.
    If necessary, the data is converted from machine independent format into
    the format needed by the platform. See |SysWriteLong()| for details.
 ** @see SysWriteLong MfReadLong
 **/

int MfWriteLong(MtxFile_t *f, const long *buf, int count)

{
    int rc;
    if (!MfIsValid(f))
	return -1;
    rc = SysWriteLong(f->File,buf,count);
    if (rc != count)
	MTX_ERROR1("%s: write failed",f->Name);
    return rc;
}
