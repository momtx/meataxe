/* ============================= C MeatAxe ==================================
   File:        $Id: matech.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Reduce a matrix to semi echelon form.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <stdlib.h>



MTX_DEFINE_FILE_INFO



static int zmkechelon(PTR matrix, int nor, int noc, int *piv, int *ispiv)

{
    PTR x, newrow;
    int i, j, rank;


    /* Initialize the table
       -------------------- */
    for (i = 0; i < noc; ++i)
    {
	piv[i] = i;
	ispiv[i] = 0;
    }

    /* Echelonize the matrix and build the pivot table in <piv>.
       Keep track of assigned pivot columns in <ispiv>.
       --------------------------------------------------------- */
    rank = 0;
    newrow = matrix;
    for (i = 0, x = matrix; i < nor && rank < noc; ++i, FfStepPtr(&x))
    {
        int newpiv;
        FEL f;

        if (rank < i)
            FfCopyRow(newrow,x);
        FfCleanRow(newrow,matrix,rank,piv);
        newpiv = FfFindPivot(newrow,&f);
	MTX_ASSERT(newpiv < noc);
        if (newpiv >= 0)
        {
	    piv[rank] = newpiv;
	    ispiv[newpiv] = 1;
	    ++rank;
            FfStepPtr(&newrow);
        }
    }

    /* Insert the non-pivot columns
       ---------------------------- */
    j = rank;
    for (i = 0; i < noc; ++i)
    {
	if (!ispiv[i])
	    piv[j++] = i;
    }
    MTX_VERIFY(j == noc);

    return rank;
}


/// @addtogroup mat
/// @{


/// Reduce to echelon form
/// This function performs a Gaussian elimination on the matrix |mat|. On
/// return, |mat| is in semi echelon form and a pivot table has been 
/// attatched to the matrix. If the rank of |mat| was smaller than the number 
/// of rows, some rows are removed during the process. This function can also 
/// be used to rebuild the pivot table after the matrix has been modified.
/// @param mat Pointer to the matrix.
/// @return Rank of @em mat, or -1 on error.

int MatEchelonize(Matrix_t *mat)
{
    int rank;
    int *newtab;
    static int *is_pivot = NULL;
    static int maxnoc = -1;

    /* Check the argument
       ------------------ */
    if (!MatIsValid(mat))
	return -1;

    /* Re-allocate the pivot table. This is not really necessary, since
       |Noc| should never change without releasing the pivot table, but
       this would be a really nasty bug....
       ----------------------------------------------------------------- */
    newtab = NREALLOC(mat->PivotTable,int,mat->Noc);
    if (newtab == NULL)
    {
	MTX_ERROR1("Cannot allocate pivot table (size %d)",mat->Noc);
	return -1;
    }
    mat->PivotTable = newtab;
    if (mat->Noc > maxnoc)
    {
	int *new_is_piv = NREALLOC(is_pivot,int,mat->Noc);
	if (new_is_piv == NULL)
	{
	    MTX_ERROR("Cannot allocate temporary table");
	    return -1;
	}
	is_pivot = new_is_piv;
	maxnoc = mat->Noc;
    }

    /* Build the pivot table
       --------------------- */
    FfSetField(mat->Field);
    FfSetNoc(mat->Noc);
    rank = zmkechelon(mat->Data,mat->Nor,mat->Noc,mat->PivotTable,is_pivot);

    /* If the rank is less than the number of rows, remove null rows
       ------------------------------------------------------------- */
    if (rank != mat->Nor)
    {
	mat->Nor = rank;
	mat->Data = (PTR) SysRealloc(mat->Data,FfCurrentRowSize * rank);
    }

    return rank;
}





/// Nullity of a matrix.
/// This function calculates the dimension of the null-space of a matrix.
/// Unlike MatNullity__() this function does not modify the matrix.
/// @param mat Pointer to the matrix.
/// @return Nullity of the matrix, or -1 on error.

long MatNullity(const Matrix_t *mat)
{
    return MatNullity__(MatDup(mat));
}



/// Nullity of a matrix.
/// This function calculates the dimension of the null-space of a matrix
/// and deletes the matrix.
/// @param mat Pointer to the matrix.
/// @return Nullity of @em mat, or -$ on error.

long MatNullity__(Matrix_t *mat)
{
    long nul;
    MatEchelonize(mat);
    nul = mat->Noc - mat->Nor;
    MatFree(mat);
    return nul;
}


/// @}
