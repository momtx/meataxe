/* ============================= C MeatAxe ==================================
   File:        $Id: setcore.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Basic set functions.
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

static const int InitialSize = 10;
static const unsigned long SetMagic = 0xEF452338;




/// @defgroup intset Sets of Integers.
/// @{

/// @class Set_t
/// @brief
/// A Set of Integers.
/// The Set_t structure represents a set of (long) integers. Internally, the set is stored
/// as a sorted list. Insert operations are relatively expensive, especially for large sets.
/// So, if you expect a lot of inserts, the BitString_t data type may be a better choice.


/// Check a set.
/// This function checks if the argument is a valid set. If the set is o.k.,
/// the function returns 1. 
/// Otherwise, an error is signaled and, if the error handler does not 
/// terminate the program, the function returns 0.
/// @param s Pointer to the set.
/// @return 1 if @a s is a valid set, 0 otherwise.

int SetIsValid(const Set_t *s)
{
    if (s == NULL)
    {
	MTX_ERROR("NULL set");
	return 0;
    }
    if (s->Magic != SetMagic || s->Size < 0 || s->BufSize < s->Size)
    {
	MTX_ERROR3("Invalid set (Magic=%d, Size=%d, BufSize=%d)",
	    (int)s->Magic,s->Size,s->BufSize);
	return 0;
    }
    if (s->Data == NULL)
    {
	MTX_ERROR("Data==NULL in set");
	return 0;
    }
    return 1;
}




/// Create a new set.
/// This function creates a new, empty set. To destroy a set, 
/// use SetFree(), @em not SysFree().
/// @return Pointer to the new set or 0 on error.

Set_t *SetAlloc()
{
    Set_t *x;

    x = ALLOC(Set_t);
    if (x == NULL)
    {
 	MTX_ERROR("Cannot allocate set");
	return NULL;
    }
    x->Size = 0;
    x->BufSize = InitialSize;
    x->Data = NALLOC(long,InitialSize);
    if (x->Data == NULL)
    {
	SysFree(x);
 	MTX_ERROR("Cannot allocate set data");
	return NULL;
    }
    x->Magic = SetMagic;
    return x;
}




/// Destroy a set.
/// This function frees an integer set. The argument must be a Set_t 
/// structure which has previously been allocated with SetAlloc(). 
/// @param x Pointer to the set.
/// @return 0 on success, -1 on error.

int SetFree(Set_t *x)
{
    if (!SetIsValid(x))
	return -1;
    SysFree(x->Data);
    memset(x,0,sizeof(Set_t));
    SysFree(x);
    return 0;
}




/// Duplicate a set.
/// @param s Pointer to the set.
/// @return Pointer to a copy of the set, or 0 on error.

Set_t *SetDup(const Set_t *s)
{
    Set_t *x;

    if (!SetIsValid(s))
	return NULL;
    x = ALLOC(Set_t);
    if (x == NULL)
    {
 	MTX_ERROR("Cannot allocate set");
	return NULL;
    }

    x->Size = s->Size;
    x->BufSize = s->Size;
    x->Data = NALLOC(long,x->BufSize);
    if (x->Data == NULL)
    {
	SysFree(x);
 	MTX_ERROR("Cannot allocate set data");
	return NULL;
    }
    memcpy(x->Data,s->Data,sizeof(x->Data[0]) * s->Size);
    x->Magic = SetMagic;
    return x;
}


/// @}
