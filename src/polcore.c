/* ============================= C MeatAxe ==================================
   File:        $Id: polcore.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Basic polynomial functions.
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

#define POLY_MAGIC 0x355A3207



/**
 ** @defgroup poly Polynomials
 ** @details
 ** The MeatAxe can work with polynomials over a finite field. A polynomial is represented
 ** by a Poly_t structure. Each polynomial carries the field order, i.e., you can work
 ** with polynomials over different fields on one program. However, this feature is 
 ** currently of little use since all standard operations only work on polynomials over the
 ** same field, and there is no easy way to identify polynomials over a field and its subfields.
 **
 ** There is a second representation of polynomials as product of factors, see FPoly_t.
 **/

/** 
 ** @addtogroup poly
 ** @{
 **/

/** @class Poly_t
 ** A Polynomial.
 ** Internally, a polynomial of degree n is represented as an array of n+1 field
 ** elements (@c Data field), where <tt>data[i]</tt> is the coefficient of x^i.
 ** The leading coefficient is always non-zero on the MeatAxe API level (it can 
 ** temporarily be zero during calculations). The null polynomial has a degree of -1.
 **/
   
/**
 ** Check a polynomial.
 ** This function checks if the argument is a pointer to a valid polynomial. If the polynomial
 ** the function returns 1. Otherwise, an error is signalled and, if the error handler does not
 ** terminate the program, the function returns 0.
 ** @param p The polynomial to check.
 ** @return 1 if @em p points to a valid polynomial, 0 otherwise.
 **/

int PolIsValid(const Poly_t *p)
{
    int deg;

    if (p == NULL)
    {
	MTX_ERROR("NULL polynomial");
	return 0;
    }
    deg = p->Degree;
    if (p->Magic != POLY_MAGIC || deg < -1 || p->Field < 2 
	|| p->Data == NULL || p->BufSize < 0)
    {
	MTX_ERROR4("Invalid polynomial (magic=%d, field=%d, deg=%d, bufsize=%d)",
	    p->Magic,p->Field,deg,p->BufSize);
	return 0;
    }

    if (deg >= 0 && p->Data[deg] == FF_ZERO)
    {
	MTX_ERROR("Invalid polynomial: leading coefficient is zero");
	return 0;
    }
    return 1;
}







/**
 ** Allocate a polynomial
 ** This function creates the polynomial p(x)=x^n over the current field.
 ** If n is negative, a zero polynomial is created. The return value is a 
 ** pointer to a newly allocated Poly_t structure. The caller is responsible 
 ** for releasing memory by calling PolFree() when the polynomial is no 
 ** longer needed.
 ** @param fl Field order.
 ** @param n Degree of the polynomial.
 ** @return Pointer to a new Poly_t structure or 0 on error.
 **/

Poly_t *PolAlloc(int fl, int n)
{
    Poly_t *x;
    int i, s;

    if (n < 0) 
	n = -1;
    if ((s = n + 1) < 1) 
	s = 1;

    FfSetField(fl);
    if ((x = ALLOC(Poly_t)) == NULL) 
    {
	MTX_ERROR("Cannot allocate polynomial");
	return NULL;
    }
    x->Magic = POLY_MAGIC;
    x->Field = fl;
    x->Degree = n;
    x->BufSize = s;
    if ((x->Data = NALLOC(FEL,s)) == NULL)
    {
	SysFree(x);
        MTX_ERROR("Cannot allocate polynomial data");
        return NULL;
    }
    for (i = 0; i < (int) s-1; ++i) 
	x->Data[i] = FF_ZERO;
    x->Data[s-1] = FF_ONE;
    return x;
}




/**
 ** Free a polynomial"
 ** This function frees a polynomial data structure and cleans up all internal data.
 ** @param x Pointer to the polynomial.
 ** @return $0$ on success, $-1$ on error.
 **/

int PolFree(Poly_t *x)

{
    if (!PolIsValid(x))
	return -1;
    SysFree(x->Data);
    memset(x,0,sizeof(Poly_t));
    SysFree(x);
    return 0;
}


/**
 ** Normalize a polynomial.
 ** This function makes sure that the leading coefficient of a polynomial is non-zero.
 **/

void Pol_Normalize(Poly_t *p)
{
    int i = p->Degree;
    while (i >= 0 && p->Data[i] == FF_ZERO) 
	--i;
    p->Degree = i;
}



/** 
 ** @}
 **/
