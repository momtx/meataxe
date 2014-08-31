/* ============================= C MeatAxe ==================================
   File:        $Id: rdcfgen.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Read generators for constituents.
   --------------------------------------------------------------------------
   (C) Copyright 1997  Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"

MTX_DEFINE_FILE_INFO

/**
 ** @addtogroup cfinfo
 ** @{
 **/

/**
 ** Load a Constituent.
 ** This function reads the generators of one constituent of a module and optionally performs
 ** some basic operations (invert, transpose) on the generators.
 ** @a flags may be any combination of the following values:
 ** - @c LAT_RG_STD: Read generators in standard basis (file XY.std.N). These files are
 **   created by pwkond. Default is to read the generators in "random" basis as they are
 **   produced by @ref prog_chop "chop".
 ** - @c LAT_RG_INVERT: Invert the generators.
 ** - @c LAT_RG_TRANSPOSE: Transpose the generators.
 **
 ** @param info Pointer to the lattice information structure.
 ** @param cf Number of constituent to load.
 ** @param flags Optional flags (see below).
 ** @return Pointer to a matrix representation of the specified constituent,
    or 0 on error.
 **/

MatRep_t *Lat_ReadCfGens(Lat_Info *info, int cf, int flags)
{
    MatRep_t *rep;
    char fn[sizeof(info->BaseName) + 20];
    int i;

    /* Check the arguments
       ------------------- */
    if (info == NULL)
    {
	MTX_ERROR1("info: %E",MTX_ERR_BADARG);
	return NULL;
    }
    if (cf < 0 || cf >= info->NCf)
    {
	MTX_ERROR1("cf: %E",MTX_ERR_BADARG);
	return NULL;
    }

    /* Make file name
       -------------- */
    if (flags & LAT_RG_STD)
	sprintf(fn,"%s%s.std.%%d",info->BaseName,Lat_CfName(info,cf));
    else
	sprintf(fn,"%s%s.%%d",info->BaseName,Lat_CfName(info,cf));

    /* Load the representation
       ----------------------- */
    rep = MrLoad(fn,info->NGen);
    if (rep == NULL)
	return NULL;

    /* Apply modifications
       ------------------- */
    for (i = 0; i < rep->NGen; ++i)
    {
	if (flags & LAT_RG_INVERT)
	{
	    Matrix_t *mat2 = MatInverse(rep->Gen[i]);
	    if (mat2 == NULL)
	    {
		MTX_ERROR1("Cannot transpose generator %d",i);
		MrFree(rep);
		return NULL;
	    }
	    MatFree(rep->Gen[i]);
	    rep->Gen[i] = mat2;
	}
	if (flags & LAT_RG_TRANSPOSE)
	{
	    Matrix_t *mat2 = MatTransposed(rep->Gen[i]);
	    if (mat2 == NULL)
	    {
		MTX_ERROR1("Cannot invert generator %d",i);
		MrFree(rep);
		return NULL;
	    }
	    MatFree(rep->Gen[i]);
	    rep->Gen[i] = mat2;
	}
    }
    return rep;
}

/**
 ** @}
 **/
