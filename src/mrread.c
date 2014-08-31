/* ============================= C MeatAxe ==================================
   File:        $Id: mrread.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Load a matrix representation.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <string.h>

MTX_DEFINE_FILE_INFO

/**
 ** @addtogroup mrep
 ** @{
 **/


/**
 ** Load a Matrix Representation.
 ** This function creates a new matrix representation and reads the generators
 ** from files. Each generator ist expected in a different file. The file name
 ** is constructed by appending ".1", ".2" etc. to @a basename or, if @a basename
 ** contains a "%d" placeholder, by replacing the "%d" with "1", "2", etc.
 ** For example, the following lines 
 ** @code
 ** m11 = MrLoad("m11",2);
 ** m11 = MrLoad("m11.%d",2);
 ** @endcode
 ** are equivalent. In both cases, two matrices are read from "m11.2" and "m11.2",
 ** repectively.
 ** @param basename Base file name for generators.
 ** @param ngen Number of generators.
 ** @return Pointer to the representation, or 0 on error.
 **/

MatRep_t *MrLoad(const char *basename, int ngen)
{
    char *fn;
    int ext_format;	    /* '%d' found in <basename> */
    MatRep_t *mr;
    int i;

    /* Make a copy of the basename and reserve extra bytes for the extension.
       ---------------------------------------------------------------------- */
    fn = SysMalloc(strlen(basename) + 10);
    if (fn == NULL)
    {
	MTX_ERROR("Cannot allocate buffer");
	return NULL;
    }

    /* Allocate the Representation
       --------------------------- */
    mr = MrAlloc(0,NULL,0);
    if (mr == NULL)
    {
	MTX_ERROR("Cannot allocate representation");
	SysFree(fn);
	return NULL;
    }

    /* Read the generators
       ------------------- */
    ext_format = strstr(basename,"%d") != NULL;
    for (i = 0; i < ngen; ++i)
    {
	Matrix_t *gen;
	if (ext_format)
	    sprintf(fn,basename,i+1);
	else
	    sprintf(fn,"%s.%d",basename,i+1);
	if ((gen = MatLoad(fn)) == NULL || MrAddGenerator(mr,gen,0) != 0)
	{
	    MTX_ERROR("Cannot load generator");
	    MrFree(mr);
	    SysFree(fn);
	    return NULL;
	}
    }

    SysFree(fn);
    return mr;
}

/**
 ** @}
 **/

