////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Matrix set core functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

   
////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

#define MS_MAGIC 0x6263659B
   

/// @defgroup matset Matrix Sets
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @class MatrixSet_t
/// @brief A Set of Matrices.
/// The MatrixSet_t structure represents a sequence of linearly independent
/// matrices over a finite field. All matrices in the set have the same number
/// of rows and columns. To assure that the matrices are linearly independent,
/// the sequence is kept in echelon form. This means, each matrix in the set
/// has a nonzero mark, and all matrices after this matrix have a zero at the
/// corresponding position.
/// Note that the individual matrices within the set are not in echelon form.
///
/// There is only one way to create a matrix set: the application first
/// allocates a MatrixSet_t structure, and then adds matrices to the 
/// set with msCleanAndAppend(). The latter function guarantees that
/// the matrix set remains linearly independent --- it will not add a 
/// matrix which is already in the span of the set.
/// A second function, msClean(), can be used to determine if a
/// matrix is in the span of a matrix set without modifying the set.
///
/// Once a matrix has been added to a matrix set, the set takes the
/// ownership of that matrix. The application must not modify of free 
/// a matrix after it has een added to a matrix set. When the matrix 
/// set is freed with msFree(), all matrices in  the set are freed, 
/// too.


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check a matrix set.
/// This function checks if the argument |set| is a pointer to a valid
/// matrix. If the matrix set is valid, the function returns 1. Otherwise, 
/// an error is signalled and, if the error handler does not terminate the 
/// program, the function returns 0.
/// @param set Pointer to the matrix set.
/// @return 1 if @a set points to a valid matrix set, 0 otherwise.

int msIsValid(const MatrixSet_t *set)
{
#ifdef PARANOID
    int i;
    MatrixSetElement_t *l;
    int field, nor, noc;
#endif

    if (set == NULL || set->Magic != MS_MAGIC)
    {
	mtxAbort(MTX_HERE,"Invalid matrix set at 0x%lx",(long) set);
	return 0;
    }
    if (set->Len < 0)
    {
	mtxAbort(MTX_HERE,"Invalid matrix set: len=%d",set->Len);
	return 0;
    }
    if (set->Len > 0 && set->List == NULL)
    {
	mtxAbort(MTX_HERE,"Invalid matrix set: list=NULL");
	return 0;
    }

#ifdef PARANOID
    for (i = 0, l = set->List; i < set->Len; ++i, ++l)
    {
	if (!matValidate(MTX_HERE, l->Matrix))
	{
	    mtxAbort(MTX_HERE,"Matrix %d is invalid",i);
	    return 0;
	}
	if (i == 0)
	{
	    field = l->Matrix->Field;
	    nor = l->Matrix->Nor;
	    noc = l->Matrix->Noc;
	}
	else if (field != l->Matrix->Field ||
	    nor != l->Matrix->Nor || noc != l->Matrix->Noc)
	{
	    mtxAbort(MTX_HERE,"Invalid matrix set: Matrix %d and %d: %s",
		i,i-1,MTX_ERR_INCOMPAT);
	    return 0;
	}
	if (l->PivRow < 0 || l->PivRow >= nor)
	{
	    mtxAbort(MTX_HERE,"Invalid pivot row %d, nor=%d",l->PivRow,nor);
	    return 0;
	}
	if (l->PivCol < 0 || l->PivCol >= noc)
	{
	    mtxAbort(MTX_HERE,"Invalid pivot column %d, noc=%d",l->PivRow,noc);
	    return 0;
	}
    }
#endif
    return 1;
}




////////////////////////////////////////////////////////////////////////////////////////////////////
/// Allocate a matrix set.
/// This function allocates a new matrix set. The set is initially empty.
/// When the set is no longer needed, if must be freed with msFree().
/// @return Pointer to the new matrix set, or 0 on error.

MatrixSet_t *msAlloc()
{
    MatrixSet_t *x;

    x = ALLOC(MatrixSet_t);
    if (x == NULL)
    {
	mtxAbort(MTX_HERE,"Cannot allocate matrix set");
	return NULL;
    }
    memset(x,0,sizeof(*x));
    x->Magic = MS_MAGIC;
    return x;
}




////////////////////////////////////////////////////////////////////////////////////////////////////
/// Free a matrix set.
/// This function frees a matrix set. All matrices in the set are freed, too.
/// @param set Pointer to the matrix set.
/// @return 0 on success, -1 on error.

int msFree(MatrixSet_t *set)
{
    int i;
    if (!msIsValid(set))
	return -1;
    for (i = 0; i < set->Len; ++i)
	matFree(set->List[i].Matrix);
    sysFree(set->List);
    memset(set,0,sizeof(*set));
    return 0;
}


/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
