/* ============================= C MeatAxe ==================================
   File:        $Id: polwrite.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Write polynomial into a file.
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



/* -----------------------------------------------------------------
   mktmp() - Allocate temporary workspace
   ----------------------------------------------------------------- */

static long tmpfl = 0;
static long tmpdeg = 0;
static PTR tmpvec = NULL;

static void mktmp(long fl, long deg)

{
    FfSetField(fl);
    if (deg > 0) FfSetNoc(deg+1);
    if (tmpfl != fl || tmpdeg < deg)
    {
	if (tmpvec != NULL) SysFree(tmpvec);
	tmpvec = FfAlloc(1);
	tmpdeg = deg;
	tmpfl = fl;
    }
}


/// @addtogroup poly
/// @{


/// Write a polynomial to a file.
/// @see PolSave()
/// @param p Pointer to the polynomial.
/// @param f File to write to.
/// @return 0 on success, -1 on error.

int PolWrite(const Poly_t *p, FILE *f)
{
    long hdr[3];

    if (!PolIsValid(p))
	return -1;
    mktmp(p->Field,p->Degree);
    hdr[0] = -2;
    hdr[1] = p->Field;
    hdr[2] = p->Degree;
    if (SysWriteLong(f,hdr,3) != 3) 
    {
	MTX_ERROR("Cannot write header");
	return -1;
    }
    if (p->Degree >= 0)
    {
	int i;
    	for (i = 0; i <= p->Degree; ++i)
	    FfInsert(tmpvec,i,p->Data[i]);
        if (FfWriteRows(f,tmpvec,1) != 1)
	{
	    MTX_ERROR("Cannot write data");
	    return -1;
	}
    }
    return 0;
}




/// Write a Polynomial to a File.
/// This function creates a file, writes a single polynomial to the file and
/// closes the file. If a f ile with the specified name already exists, it's
/// contents are destroyed.
/// @see PolWrite
/// @param pol Polynomial to write.
/// @param fn File name.
/// @return 0 on success, -1 on error.

int PolSave(const Poly_t *pol, const char *fn)

{
    FILE *f;
    int result;

    if (!PolIsValid(pol))
	return -1;
    if ((f = SysFopen(fn,FM_CREATE)) == NULL)
    {
	MTX_ERROR1("Cannot open %s",fn);
	return -1;
    }
    result = PolWrite(pol,f);
    fclose(f);
    if (result != 0)
    {
	MTX_ERROR1("Cannot write polynomial to %s",fn);
	return -1;
    }
    return result;
}



/// @}
