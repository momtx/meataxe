////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Matrix representations, core functions
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
///   Local data
   
MTX_DEFINE_FILE_INFO

#define MR_MAGIC 0x1bb50442

   
/// @defgroup mrep Matrix Representations
/// @{
/// @class MatRep_t
/// A matrix representation.
/// This data structure is basically a collection of N matrices, representing the generators
/// of the algebra. A matrix representation is created with MrAlloc(). If the generators are
/// available, they can be passed to MrAlloc(). Generators can also be added to a matrix
/// representation with MrAddGenerator(). In both cases the user can choose if a copy of the
/// matrices or only a reference to the matrices is stored in the representation. In either
/// case, deleting the representation with MrFree() will also delete the generators.
/// For example, after the following code has been executed:
/// @code
/// Matrix_t *mat = MatLoad("matrix");
/// Rep = MrAlloc(1,&mat,0);
/// MrFree(Rep);                               // *** invalidates mat
/// @endcode
/// the @c mat pointer is no longer valid, because the matrix was deleted by MrFree().
/// However, after
/// @code
/// Matrix_t *mat = MatLoad("matrix");
/// Rep = MrAlloc(1,&mat,MR_COPY_GENERATORS);
/// MrFree(Rep);                               // *** mat remains valid
/// @endcode
/// @c mat is still valid because a copy of the matrix was created.

////////////////////////////////////////////////////////////////////////////////////////////////////

static int GensAreValid(int ngen, Matrix_t **gen)
{
    int i;

    if (ngen < 0)
    {
	MTX_ERROR1("ngen: %E",MTX_ERR_BADARG);
	return 0;
    }
    if (ngen > 0 && gen == NULL)
    {
	MTX_ERROR1("gen == NULL: %E",MTX_ERR_BADARG);
	return 0;
    }
    for (i = 0; i < ngen; ++i)
    {
	if (!MatIsValid(gen[i]))
	{
	    MTX_ERROR1("gen[%d] invalid",i);
	    return 0;
	}
	if (gen[i]->Nor != gen[i]->Noc)
	{
	    MTX_ERROR2("gen[%i]: %E",i,MTX_ERR_NOTSQUARE);
	    return 0;
	}
	if (i != 0)
	{
	    if (gen[i]->Field != gen[0]->Field || gen[i]->Nor != gen[0]->Nor)
	    {
		MTX_ERROR2("gen[0] and gen[%d]: %E",i,MTX_ERR_INCOMPAT);
		return 0;
	    }
	}
    }
    return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check a matrix representation.
/// This function checks if the argument is a pointer to a valid
/// matrix representation. If the representation is o.k., the function 
/// returns 1. Otherwise, an error is signaled and, if the error handler 
/// does not terminate the program, the function returns 0.
/// @param rep Pointer to the matrix representation.
/// @return 1 if @a rep points to a valid matrix representation, 0 otherwise.

int MrIsValid(const MatRep_t *rep)
{
    if (rep == NULL)
    {
	MTX_ERROR("NULL representation");
	return 0;
    }
    if (rep->Magic != MR_MAGIC)
    {
	MTX_ERROR1("Invalid matrix representation (magic=%d)",(int)rep->Magic);
	return 0;
    }
    if (!GensAreValid(rep->NGen,rep->Gen))
    {
	MTX_ERROR("Invalid generators");
	return 0;
    }
    return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Create a matrix representation.
/// This function creates a new matrix representation for a given set of generators.
/// The matrices in @a gen must all be  square, over the same field, and with the same dimensions.
/// @a flags may be zero or the special value MR_COPY_GENERATORS. In the latter case,
/// a local copy of the generators is made, and the matrices in @a gen can
/// be safely destroyed. If @a flags is 0, only references to the matrices 
/// are stored in the MatRep_t structure. Consequently, the application 
/// must not modify or destroy the matrices after calling MrAlloc(). They
/// will be destroyed automatically when MrFree() is called to destroy the
/// representation.
/// @param ngen Number of generators in @a gen.
/// @param gen List of generators.
/// @param flags Optional flags (see function description).
/// @return Pointer to the new matrix representation or 0 on error.

// TODO: add a MR_FOREIGN_GENERATORS to prevent MrFree() from destroying the generators

MatRep_t *MrAlloc(int ngen, Matrix_t **gen, int flags)
{
    MatRep_t *rep;
    int i;

    if (!GensAreValid(ngen,gen))
    {
	MTX_ERROR1("%E",MTX_ERR_BADARG);
	return NULL;
    }

    // Allocate a new MatRep_t structure
    rep = ALLOC(MatRep_t);
    if (rep == NULL)
    {
	MTX_ERROR("Cannot allocate MatRep_t structure");
	return NULL;
    }
    memset(rep,0,sizeof(MatRep_t));
    rep->Gen = NALLOC(Matrix_t *,ngen);
    if (rep->Gen == NULL)
    {
	MTX_ERROR("Cannot allocate generator list");
	SysFree(rep);
	return NULL;
    }

    // Copy generators
    rep->NGen = ngen;
    for (i = 0; i < ngen; ++i)
    {
	if (flags & MR_COPY_GENERATORS)
	{
	    rep->Gen[i] = MatDup(gen[i]);
	    if (rep->Gen[i] == NULL)
	    {
		MTX_ERROR("Cannot copy generator");
		while (--i >= 0) {
		    MatFree(rep->Gen[i]);
		}
		SysFree(rep->Gen);
		SysFree(rep);
		return NULL;
	    }
	}
	else
	    rep->Gen[i] = gen[i];
    }

    rep->Magic = MR_MAGIC;
    return rep;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Delete a matrix representation.
/// This function frees a matrix representation which has beed created by 
/// MrAlloc(). This implies freeing the internal data buffers as well as 
/// the MatRep_t structure itself. Note: Even if the representation 
/// was created without MR_COPY_GENERATORS, the matrices that were passed 
/// to MrAlloc() are now destroyed. The same applies to matrices added to 
/// the representation with MrAddGenerator().
/// @param rep Pointer to the matrix representation.
/// @return 0 on success, -1 on error.

int MrFree(MatRep_t *rep)
{
    int i;
    if (!MrIsValid(rep))
    {
	MTX_ERROR1("%E",MTX_ERR_BADARG);
	return -1;
    }
    for (i = 0; i < rep->NGen; ++i)
        MatFree(rep->Gen[i]);
    memset(rep->Gen,0,sizeof(Matrix_t *) * rep->NGen);
    SysFree(rep->Gen);
    memset(rep,0,sizeof(MatRep_t));
    SysFree(rep);
    return 0;
}

/// @}
