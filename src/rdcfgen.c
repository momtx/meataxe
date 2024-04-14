////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Read generators for constituents.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"


/// @addtogroup cfinfo
/// @{

/// Load a Constituent.
/// This function reads the generators of one constituent of a module and optionally performs
/// some basic operations (invert, transpose) on the generators.
/// @p flags may be any combination of the following values:
/// - @c LAT_RG_STD: Read generators in standard basis (file XY.std.N). These files are
///   created by pwkond. Default is to read the generators in "random" basis as they are
///   produced by @ref prog_chop "chop".
/// - @c LAT_RG_INVERT: Invert the generators.
/// - @c LAT_RG_TRANSPOSE: Transpose the generators.
///
/// @param info Pointer to the lattice information structure.
/// @param cf Number of constituent to load.
/// @param flags Optional flags (see below).
/// @return Pointer to a matrix representation of the specified constituent,
/// or 0 on error.

MatRep_t *latReadCfGens(LatInfo_t *info, int cf, int flags)
{
    if (info == NULL)
	mtxAbort(MTX_HERE,"info: %s",MTX_ERR_BADARG);
    if (cf < 0 || cf >= info->nCf)
	mtxAbort(MTX_HERE,"cf: %s",MTX_ERR_BADARG);

    // Load the representation
    MatRep_t *rep;
    {
       StrBuffer_t* fn = sbAlloc(100);
       sbAppend(fn, info->baseName);
       sbAppend(fn, latCfName(info,cf));
       if (flags & LAT_RG_STD) sbAppend(fn,".std");
       sbAppend(fn, ".%d");
       rep = mrLoad(sbData(fn),info->NGen);
       sbFree(fn);
    }

    // Apply modifications
    for (int i = 0; i < rep->NGen; ++i)
    {
	if (flags & LAT_RG_INVERT)
	{
	    Matrix_t *mat2 = matInverse(rep->Gen[i]);
	    matFree(rep->Gen[i]);
	    rep->Gen[i] = mat2;
	}
	if (flags & LAT_RG_TRANSPOSE)
	{
	    Matrix_t *mat2 = matTransposed(rep->Gen[i]);
	    matFree(rep->Gen[i]);
	    rep->Gen[i] = mat2;
	}
    }
    return rep;
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
