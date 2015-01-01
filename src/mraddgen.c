////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Add a generator to a matrix representation
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"

MTX_DEFINE_FILE_INFO

/// @addtogroup mrep
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Add a Matrix to a Representation.
/// This function adds a generator, i.e., a matrix, to a matrix representation. 
/// The matrix must be square. If there are already generators in the 
/// represenatation, the new matrix must be over the same field and have 
/// the same number of rows. @a flags may be zero or the special value 
/// MR_COPY_GENERATORS. In the latter case, a local copy of the generator 
/// is made, and @a gen can be safely destroyed. Otherwise, only a reference to 
/// the matrix is stored in the MatRep_t structure. Consequently, the application
/// must not modify or destroy the matrix after calling %MrAddGenerator(). 
/// It will be destroyed by MrFree() together with the representation.
/// @param rep Matrix representation.
/// @param gen Matrix to add.
/// @param flags Optional flags (see below).
/// @return 0 on success, -1 on error.

int MrAddGenerator(MatRep_t *rep, Matrix_t *gen, int flags)
{
    Matrix_t **new_gen;

    /* Check arguments
       --------------- */
    if (!MrIsValid(rep))
    {
	MTX_ERROR1("rep: %E",MTX_ERR_BADARG);
	return -1;
    }
    if (gen->Nor != gen->Noc)
    {
	MTX_ERROR1("gen: %E",MTX_ERR_NOTSQUARE);
	return -1;
    }
    if (rep->NGen > 0)
    {
	if (gen->Field != rep->Gen[0]->Field || gen->Nor != rep->Gen[0]->Nor)
	{
	    MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
	    return -1;
	}
    }

    /* Extend the matrix list
       ---------------------- */
    new_gen = NREALLOC(rep->Gen,Matrix_t *,rep->NGen + 1);
    if (new_gen == NULL)
    {
	MTX_ERROR("Cannot extend matrix list");
	return -1;
    }

    /* Make a copy of the generator, if requested
       ------------------------------------------ */
    if (flags & MR_COPY_GENERATORS)
    {
	gen = MatDup(gen);
        if (gen == NULL)
	{
	    MTX_ERROR("Cannot copy generator");
	    return -1;
	}
    }

    new_gen[rep->NGen] = gen;
    rep->Gen = new_gen;
    ++rep->NGen;

    return 0;
}


/// @}

