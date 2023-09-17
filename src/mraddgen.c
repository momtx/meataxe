////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Add a generator to a matrix representation
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"


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
/// must not modify or destroy the matrix after calling %mrAddGenerator(). 
/// It will be destroyed by mrFree() together with the representation.
/// @param rep Matrix representation.
/// @param gen Matrix to add.
/// @param flags Optional flags (see below).
/// @return 0 on success, -1 on error.

int mrAddGenerator(MatRep_t *rep, Matrix_t *gen, int flags)
{
    Matrix_t **new_gen;

    /* Check arguments
       --------------- */
    if (!mrIsValid(rep))
    {
	mtxAbort(MTX_HERE,"rep: %s",MTX_ERR_BADARG);
    }
    if (gen->Nor != gen->Noc)
    {
	mtxAbort(MTX_HERE,"gen: %s",MTX_ERR_NOTSQUARE);
    }
    if (rep->NGen > 0)
    {
	if (gen->Field != rep->Gen[0]->Field || gen->Nor != rep->Gen[0]->Nor)
	{
	    mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
	}
    }

    /* Extend the matrix list
       ---------------------- */
    new_gen = NREALLOC(rep->Gen,Matrix_t *,rep->NGen + 1);
    if (new_gen == NULL)
    {
	mtxAbort(MTX_HERE,"Cannot extend matrix list");
    }

    /* Make a copy of the generator, if requested
       ------------------------------------------ */
    if (flags & MR_COPY_GENERATORS)
    {
	gen = matDup(gen);
        if (gen == NULL)
	    mtxAbort(MTX_HERE,"Cannot copy generator");
    }

    new_gen[rep->NGen] = gen;
    rep->Gen = new_gen;
    ++rep->NGen;
    return 0;
}


/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
