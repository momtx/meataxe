/* ============================= C MeatAxe ==================================
   File:        $Id: matins.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Insert a matrix into a polynomial
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"

MTX_DEFINE_FILE_INFO


/**
 ** @addtogroup mat
 ** @{
 **/

/**
 ** Insert a matrix into a polynomial
 ** Given a square matrix A and a polynomial p over the same field, this functions
 ** calculates p(A). Unlike MatInsert() this function is destructive. The result
 ** is stored in the original matrix and the old value is lost.
 ** @param mat Pointer to the matrix.
 ** @param pol Pointer to the polynomial.
 ** @return The function returns @em mat, or 0 on error.
 **/

Matrix_t *MatInsert_(Matrix_t *mat, const Poly_t *pol)
{
    Matrix_t *x = NULL;
    int i;
    int nor;
    int l;
    PTR v;
    FEL f;

    /* Check the arguments
       ------------------- */
    if (!MatIsValid(mat))
	return NULL;
    if (!PolIsValid(pol))
	return NULL;
    if ((nor = mat->Nor) != mat->Noc) 
    {
	MTX_ERROR1("%E",MTX_ERR_NOTSQUARE);
	return NULL;
    }
    if (mat->Field != pol->Field) 
    {
	MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
	return NULL;
    }

    FfSetField(mat->Field);
    FfSetNoc(nor);

    /* Special case: p(x) = 0
       ---------------------- */
    if (pol->Degree == -1)
    {
	for (l = 0, v = mat->Data; l < nor; FfStepPtr(&v), ++l)
	    FfMulRow(v,FF_ZERO);
	return mat;
    }

    /* Special case: deg(p) = 0
       ------------------------ */
    if (pol->Degree == 0)
    {
	for (l = 0, v = mat->Data; l < nor; FfStepPtr(&v), ++l)
	{
	    FfMulRow(v,FF_ZERO);
	    FfInsert(v,l,pol->Data[0]);
	}
	return mat;
    }

    /* Evaluate p(A)
       ------------- */
    if (pol->Degree > 1) 
	x = MatDup(mat);
    if ((f = pol->Data[pol->Degree]) != FF_ONE)
    {
	for (l = nor, v = mat->Data; l > 0; --l, FfStepPtr(&v))
	    FfMulRow(v,f);
    }
    for (i = pol->Degree - 1; i >= 0; --i)
    {
        if ((f = pol->Data[i]) != FF_ZERO)
        {
	    for (l = 0, v = mat->Data; l < nor; ++l, FfStepPtr(&v))
		FfInsert(v,l,FfAdd(FfExtract(v,l),f));
        }
	if (i > 0) 
	    MatMul(mat,x);
    }
    if (pol->Degree > 1) MatFree(x);
    return mat;
}



/**
 ** Insert a matrix into a polynomial
 ** Given a square matrix A and a polynomial p over the same field, this functions
 ** calculates p(A). Unlike MatInsert_() this function returns a new matrix and does
 ** not modify the original matrix.
 ** @param mat Pointer to the matrix.
 ** @param pol Pointer to the polynomial.
 ** @return @em pol(@em mat), or 0 on error.
 **/

Matrix_t *MatInsert(const Matrix_t *mat, const Poly_t *pol)
{
    Matrix_t *x;
    int i;
    int nor;
    int l;
    PTR v;
    FEL f;

    if (!MatIsValid(mat))
	return NULL;
    if (!PolIsValid(pol))
	return NULL;
    if ((nor = mat->Nor) != mat->Noc) 
    {
	MTX_ERROR1("%E",MTX_ERR_NOTSQUARE);
	return NULL;
    }
    if (mat->Field != pol->Field) 
    {
	MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
	return NULL;
    }

    /* Special case: p = 0
       ------------------- */
    if (pol->Degree == -1)
	return MatAlloc(mat->Field,nor,nor);

    /* Special case: deg(p) = 0
       ------------------------ */
    if (pol->Degree == 0)
    {
	x = MatAlloc(mat->Field,nor,nor);
	for (l = 0, v = x->Data; l < nor; ++l, FfStepPtr(&v))
	    FfInsert(v,l,pol->Data[0]);
	return x;
    }

    /* Evaluate p(A)
       ------------- */
    x = MatDup(mat);
    if ((f = pol->Data[pol->Degree]) != FF_ONE)
    {
	for (l = nor, v = x->Data; l > 0; --l, FfStepPtr(&v))
	    FfMulRow(v,f);
    }
    for (i = pol->Degree - 1; i >= 0; --i)
    {
        if ((f = pol->Data[i]) != FF_ZERO)
        {
	    for (l = 0, v = x->Data; l < nor; ++l, FfStepPtr(&v))
		FfInsert(v,l,FfAdd(FfExtract(v,l),f));
        }
	if (i > 0) 
	    MatMul(x,mat);
    }
    return x;
}

/**
 ** @}
 **/
