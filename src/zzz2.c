/* ============================= C MeatAxe ==================================
   File:        $Id: zzz2.c,v 1.2 2007-09-09 21:38:11 mringe Exp $
   Comment:     Basic row operations.
   --------------------------------------------------------------------------
   (C) Copyright 1997 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <string.h>
#include <stdlib.h>

#define LPR (FfCurrentRowSize/sizeof(long))

/* --------------------------------------------------------------------------
   Lokal data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO


/**
 ** @addtosection ff
 ** @{
 **/



/**
 ** Characteristic of the current field.
 ** Like FfOrder, this variable may be used anywhere, but it must not be modified directly.
 **/

int FfChar = 0;


/**
 ** Field order.
 ** %FfOrder may be used in expressiond but must never modified directly. To change the
 ** current field, use FfSetField().
 **/

int FfOrder = -1;


/**
 ** Field generator.
 ** This variable contains a genrator for the multiplicative group of the current field.
 **/

FEL FfGen = 0;



/**
 ** Current row size.
 ** Used by all low-level row operations. %FfNoc is updated automatically when the row size
 ** is changed with FfSetNoc().
 **/

int FfNoc = 0;



/**
 ** Allocate memory and initialize
 ** This function allocates a block of memory for a vector (if @em nrows is 1)
 ** or a matrix. Memory is initialized to zero. Field order and row size must 
 ** have been set with FfSetField() and FfSetNoc(), respectively. 
 ** @em nrows may be zero zero, i which case the functin returns a memory block of
 ** size zero which must be freed using FfFree().
 ** @param nrows Number of rows.
 ** @return Pointer to the memory block or NULL on error.
 **/

PTR FfAlloc(int nrows)
{
    PTR p;
    char *q;
    register long i;
    size_t req = (size_t) FfCurrentRowSize * (size_t) nrows;

    if ((p = (PTR) SysMalloc(req)) == NULL)
    {
	MTX_ERROR3("%E: Cannot allocate %d rows (%l bytes)",MTX_ERR_NOMEM,
	    nrows,(long)req);
	return NULL;
    }

    /* Initialize */
    q = (char *) p;
    for (i = nrows; i > 0; --i)
    {
	FfMulRow((PTR) q,FF_ZERO);
	q += FfCurrentRowSize;
    }
    return p;
}


/**
 ** Free memory.
 ** This function frees a block of memory that has been allocated with
 ** FfAlloc() before. If the argument is 0, %FfFree() does nothing.
 ** @param x Pointer to memory.
 **/

void FfFree(PTR x)
{
    if (x != NULL)
	free(x);
}



/**
 ** Copy a row.
    This function copies the contents of one row to another row.
    As with all row operations, the row length must have been set before 
    with |FfSetNoc()|.
 ** @param dest Pointer to the destination.
 ** @param src Pointer to the source.
 **/

void FfCopyRow(PTR dest, PTR src)

{
    register long *s = (long *) src, *d = (long *) dest, i;
    for (i = LPR; i > 0; --i)
   	 *d++ = *s++;
}



/**
 ** Exchange two rows
 ** This function exchanges the contents of two rows. As with all row 
 ** operations, the row length must have been set before with |FfSetNoc()|.
 ** @param dest Pointer to the first row
 ** @param src Pointer to the second row
 **/

void FfSwapRows(PTR dest, PTR src)

{	register long *p1 = (long *) src, *p2 = (long *) dest, i, l;

	for (i = LPR; i > 0; --i)
	{	l = *p1;
		*p1++ = *p2;
		*p2++ = l;
	}
}



/**
 ** Get row pointer.
 ** This function returns a pointer to the |nor|-th row of a matrix, assuming
 ** the current row size. |base| must be a pointer to the beginning of a row, 
 ** but this need not be the first row of the matrix. For example, the 
 ** following code steps through the odd rows of a matrix:
 ** @code
 ** PTR r = matrix;
 ** int i;
 ** for (i = 1; i < nrows; i += 2)
 ** {
 **    r = FfGetPtr(r,2);
 **    ...
 ** }
 ** @endcode
 ** Note: The function does not check if the resulting pointer
 ** is still inside the matrix. 
 ** @see FfStepPtr()
 ** @param base Pointer to the first row of the matrix.
 ** @param row Row index.
 **/

PTR FfGetPtr(PTR base, int row)
{
#ifdef PARANOID
    if (row < 0)
	MTX_ERROR1("row=%d < 0",row);
#endif
    return (PTR) ((char *)base + FfCurrentRowSize * row);
}


/**
 ** Advance to next row.
 ** This function increments a pointer by 1 row. The row size must have been
 ** set before with FfSetNoc(). FfStepPtr(&x) is equivalent to
 ** x = FfGetPtr(x,1). It is typically used to step through the rows of 
 ** a matrix. Here is an example with a 100 by 40 matrix:
 ** @code
 ** PTR r, mat;
 ** int i;
 ** FfSetNoc(40);
 ** mat = FfAlloc(100);
 ** for (i = 1, r = mat; i < 100; ++i, FfStepPtr(&r))
 **     ...
 ** @endcode
 ** @see FfGetPtr()
 ** @param x Pointer to the row pointer.
 ** @return Always 0.
 **/

int FfStepPtr(PTR *x)

{
    *(char **)x += FfCurrentRowSize;
    return 0;
}



/**
 ** @addtosection ff
 ** @}
 **/

