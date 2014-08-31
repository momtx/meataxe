/* ============================= C MeatAxe ==================================
   File:        $Id: temap.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Map undfer tensor product.
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"


/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO 


/**
 ** @addtogroup tp
 ** @{
 **/


/**
 ** Map Under Tensor Product.
 ** This function applies the tensor product of two matrices to one or more 
 ** vectors. The same calculation could be done with MatMul() and
 ** MatTensor(), but this function is usually faster and uses less memory,
 ** because it does not calculate the full tensor product of a⊗b.
 ** @see VectorToMatrix() MatrixToVector()
 ** @param vec Vectors to map.
 ** @param a Left matrix.
 ** @param b Right matrix.
 ** @return Image of @a vec under @a a⊗@a b, or 0 on error.
 **/

Matrix_t *TensorMap(Matrix_t *vec, const Matrix_t *a, const Matrix_t *b)
{
    Matrix_t *result;
    int i;

    /* Check the arguments.
       -------------------- */
    if (!MatIsValid(vec))
    {
	MTX_ERROR1("vec: %E",MTX_ERR_BADARG);
	return NULL;
    }
    if (!MatIsValid(a))
    {
	MTX_ERROR1("a: %E",MTX_ERR_BADARG);
	return NULL;
    }
    if (!MatIsValid(b))
    {
	MTX_ERROR1("b: %E",MTX_ERR_BADARG);
	return NULL;
    }
    if (a->Field != b->Field || b->Field != vec->Field ||
	vec->Noc != a->Nor * b->Nor)
    {
	MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
	return NULL;
    }

    /* Calculate the result.
       --------------------- */
    result = MatAlloc(vec->Field,vec->Nor,a->Noc*b->Noc);
    if (result == NULL)
	return NULL;
    for (i = 0; i < vec->Nor; ++i)
    {
	Matrix_t *tmp = MatTransposed(a);
	Matrix_t *v = VectorToMatrix(vec,i,b->Nor);
	if (v == NULL)
	{
	    MTX_ERROR("Conversion failed");
	    break;
	}
	MatMul(tmp,v);
	MatFree(v);
	MatMul(tmp,b);
	if (MatrixToVector(tmp,result,i))
	    MTX_ERROR("Conversion failed");
	MatFree(tmp);
    }
    return result;
}


/**
 ** @}
 **/
