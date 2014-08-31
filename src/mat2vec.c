/* ============================= C MeatAxe ==================================
   File:        $Id: mat2vec.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Convert matrix to vector.
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
 ** Convert matrix to vector.
 ** This function converts a matrix into a row vector by concatenating
 ** the rows of the matrix. If @a mat is a r by c matrix, the resulting
 ** vector has rc entries. Instead of allocating a new buffer for the
 ** result, %MatToVec() expects a pointer to a matrix, @a vecs, and puts 
 ** the vector into the @a n-th row of this matrix. Of course, @a vecs must
 ** be over the smae field as @a mat, have rc columns and at least n+1 rows.
 ** @see VectorToMatrix
 ** @param mat Matrix to convert.
 ** @param vecs Destination of the vector.
 ** @param n Row number where the vector is stored.
 ** @return 0 on success, -1 on error.
 **/

int MatrixToVector(const Matrix_t *mat, Matrix_t *vecs, int n)
{
    int i;

    /* Check arguments.
       ---------------- */
    if (!MatIsValid(mat))
    {
	MTX_ERROR1("mat: %E",MTX_ERR_BADARG);
	return -1;
    }
    if (!MatIsValid(vecs))
    {
	MTX_ERROR1("vecs: %E",MTX_ERR_BADARG);
	return -1;
    }
    if (   mat->Nor * mat->Noc != vecs->Noc
	|| mat->Field != vecs->Field)
    {
	MTX_ERROR1("mat and vecs: %E",MTX_ERR_INCOMPAT);
	return -1;
    }
    if (n < 0 || n >= vecs->Nor)
    {
	MTX_ERROR3("n=%d (nor=%d): %E",n,vecs->Nor,MTX_ERR_BADARG);
	return -1;
    }

    /* Convert the matrix.
       ------------------- */
    for (i = 0; i < mat->Nor; ++i)
    {
	if (MatCopyRegion(vecs,n,i*mat->Noc, mat,i,0,1,mat->Noc) != 0)
	{
	    MTX_ERROR("Copying failed");
	    return -1;
	}
    }

    return 0;
}


/**
 ** @}
 **/
