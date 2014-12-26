/* ============================= C MeatAxe ==================================
   File:        $Id: permcore.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Basic permutation functions.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"
#include <string.h>

   
/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO

#define PERM_MAGIC 0x30f8326b

   
/// @defgroup perm Permutations
/// @details
/// In the MeatAxe, a permutation of degree n operates on {0,1,...,n-1} and is represented
/// by a Perm_t structure.
/// However, in the textual representation produced by PermPrint() or by the 
/// @ref prog_zpr "zpr" program, the points are numbered from 1...n.
///
/// Only permutations of equal degree can be multiplied. This can be confusing
/// because the textual representation produced by PermPrint() does not include the
/// degree, and fixed points are always suppressed. For example "(3 4)(5 6 7)" could
/// be interpreted as a permutation of degree 8 or any higher degree. All these permutations
/// are in some natural way equal to each other, but they are different and incompatible
/// in the MeatAxe.
/// 
/// Permutations are usually created with PermAlloc() or read from a file with PermRead().
/// When a permutation is no longer used, the application must release the associated memory
/// by calling PermFree().
/// @{

/// @class Perm_t
/// @details
/// Internally, a permutation is represented as an array of long integers containing the
/// images of 0,1,...,n-1. Theoretically, the maximum degree is the largest number that
/// can be stored in a long integer. However, the MeatAxe file format is restricted to 
/// 32-bit integers, and thus the largest possible degree is 2^32.


////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Check a permutation.
/// This function checks if the argument is a pointer to a valid permutation.
/// If the permutation is o.k., the function returns 1. 
/// Otherwise, an error is signalled and, if the error handler does not 
/// terminate the program, the function returns 0.
/// @return 1 if @a p is a valid permutation, 0 otherwise.
///

int PermIsValid(const Perm_t *p)
{
    int i, deg;
    long *x;

    if (p == NULL)
    {
	MTX_ERROR("NULL permutation");
	return 0;
    }
    deg = p->Degree;
    if (p->Magic != PERM_MAGIC || deg < 0 || p->Data == NULL)
    {
	MTX_ERROR2("Invalid permutation (magic=%d, deg=%d)",
	    p->Magic,deg);
	return 0;
    }
    for (i = p->Degree, x = p->Data; i > 0; --i, ++x)
    {
	if (*x < 0 || *x >= deg)
	{
	    MTX_ERROR2("Invalid value %d in permutation (deg = %d)",
		(int) *x,deg);
	    return 0;
	}
    }

    return 1;

}


/// Allocate a permutation
/// This function creates a permutation of the specified degree.
/// The new permutation is initialized to the identity.
/// @param deg Degree.
/// @return Pointer to the new permutation, or 0 on error.

Perm_t *PermAlloc(int deg)
{
    Perm_t *p;
    int i;

    if (deg < 0)
    {
	MTX_ERROR2("deg=%d: %E",deg,MTX_ERR_BADARG);
	return NULL;
    }

    p = ALLOC(Perm_t);
    if (p == NULL)
    {
	MTX_ERROR("Cannot allocate Perm_t structure");
	return NULL;
    }
    p->Magic = PERM_MAGIC;
    p->Degree = deg;
    p->Data = NALLOC(long,deg);
    if (p->Data == NULL)
    {
	SysFree(p);
	MTX_ERROR("Cannot allocate permutation data");
	return NULL;
    }
    for (i = 0; i < deg; ++i)
	p->Data[i] = i;
    return p;
}



/// Free a permutation.
/// This function deletes a permutation and returns the memory to the 
/// system. Do not use SysFree() on permutations because this would only 
/// free the Perm_t structure but not the data buffer.
/// @param p Pointer to the permutation.
/// @return 0 on success, -1 on error.

int PermFree(Perm_t *p)
{
    /* Check arguments
       --------------- */
    if (!PermIsValid(p))
	return -1;
  
    SysFree(p->Data);
    memset(p,0,sizeof(Perm_t));
    SysFree(p);
    return 0;
}


/// @}

