/* ============================= C MeatAxe ==================================
   File:        $Id: imatcore.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Basic integer matrix functions.
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"
#include <string.h>

   
/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO

#define IMAT_MAGIC 0x396AA2F2

   
/**
 ** @defgroup imat Integer Matrices
 ** @{
 **/
   

/**
 ** @class IntMatrix_t
 ** An Integer Matrix.
 ** The IntMatrix_t structure represents a matrix with integer entries.
 ** Both @c Nor and @c Noc may be zero. In this case, @c Data ist still a valid pointer,
 ** but the memory block it points to has size zero.
 **/


/**
 ** Check an Integer Matrix.
 ** This function checks if the argument is a pointer to a valid
 ** integer matrix. If the matrix is o.k., the function returns 1. Otherwise, 
 ** an error is signalled and, if the error handler does not terminate the 
 ** program, the function returns 0.
 ** @param mat Pointer to the matrix.
 ** @return 1 if @a mat points to a valid matrix, 0 otherwise.
 **/

int ImatIsValid(const IntMatrix_t *mat)
{
    if (mat == NULL)
    {
	MTX_ERROR("NULL matrix");
	return 0;
    }
    if (mat->Magic != IMAT_MAGIC || mat->Nor < 0 || mat->Noc < 0)
    {
	MTX_ERROR2("Invalid matrix (nor=%d, noc=%d)",mat->Nor,mat->Noc);
	return 0;
    }
    return 1;

}




/**
 ** Create a New Integer Matrix.
 ** This function creates a new integer matrix with |nor| rows and |noc|
 ** columns.
 ** To destroy an integer matrix, use ImatFree(), not SysFree().
 ** @param nor Number of rows.
 ** @param noc Number of columns.
 ** @return Pointer to the new matrix or 0 on error.
 **/

IntMatrix_t *ImatAlloc(int nor, int noc)
{
    IntMatrix_t *m;

    MTX_VERIFY(nor >= 0);
    MTX_VERIFY(noc >= 0);

    /* Allocate a new Matrix_t structure
       --------------------------------- */
    m = ALLOC(IntMatrix_t);
    if (m == NULL)
    {
	MTX_ERROR("Cannot allocate IntMatrix_t structure");
	return NULL;
    }

    /* Initialize the data structure
       ----------------------------- */
    m->Magic = IMAT_MAGIC;
    m->Nor = nor;
    m->Noc = noc;
    m->Data = NALLOC(long,nor * noc);
    if (m->Data == NULL)
    {
	SysFree(m);
	MTX_ERROR("Cannot allocate matrix data");
	return NULL;
    }
    return m;
}








/**
 ** Delete an Integer Matrix.
 ** This function frees a matrix which has beed created by ImatAlloc(). This
 ** implies freeing the internal data buffers as well as the IntMatrix_t 
 ** structure itself.
 ** @param mat Pointer to the matrix.
 ** @return 0 on success, -1 on error.
 **/

int ImatFree(IntMatrix_t *mat)
{
    if (!ImatIsValid(mat))
	return -1;
    if (mat->Data != NULL)
	SysFree(mat->Data);
    memset(mat,0,sizeof(IntMatrix_t));
    SysFree(mat);
    return 0;
}


/**
 ** @}
 **/

