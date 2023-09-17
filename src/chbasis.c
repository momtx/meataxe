////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Change basis.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"



/// @addtogroup mrep 
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Change basis.
/// This function performs a basis change on a matrix representation.
/// The transformed generators are given by g'<sub>i</sub> = T g<sub>i</sub> T<sup>-1</sup>.
///
/// @param rep Matrix representation.
/// @param trans Transformation matrix mapping the old basis to the new basis. In other words,
///    the rows of the matrix are the new basis vectors

void mrChangeBasis(MatRep_t *rep, const Matrix_t *trans)
{
    Matrix_t *bi;
    int i;

    // Check arguments
    if (!mrIsValid(rep))
	mtxAbort(MTX_HERE,"rep: %s",MTX_ERR_BADARG);
    matValidate(MTX_HERE, trans);

    if (rep->NGen <= 0)
	return;
    if (trans->Field != rep->Gen[0]->Field || 
	trans->Nor != rep->Gen[0]->Nor ||
	trans->Noc != rep->Gen[0]->Noc)
    {
	mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
    }

    // Basis transformation
    if ((bi = matInverse(trans)) == NULL) 
    {
	mtxAbort(MTX_HERE,"Basis transformation is singular");
    }
    for (i = 0; i < rep->NGen; ++i)
    {
	Matrix_t *tmp = matDup(trans);
	matMul(tmp,rep->Gen[i]);
	matMul(tmp,bi);
        matFree(rep->Gen[i]);
	rep->Gen[i] = tmp;
    }
    matFree(bi);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Performs a basis change and returns a new matrix representation with the transformed
/// generators. See also @ref mrChangeBasis.

MatRep_t* mrChangeBasis2(const MatRep_t *rep, const Matrix_t *trans)
{
   if (!mrIsValid(rep))
      mtxAbort(MTX_HERE,"rep: %s",MTX_ERR_BADARG);
   matValidate(MTX_HERE, trans);
   MatRep_t *result = mrAlloc(0, NULL, 0);
   if (rep->NGen <= 0)
      return result;
   Matrix_t* g0 = rep->Gen[0];
   if (trans->Field != g0->Field || trans->Nor != g0->Nor || trans->Noc != g0->Noc)
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
   Matrix_t* bi = matInverse(trans);
   for (int i = 0; i < rep->NGen; ++i)
   {
      Matrix_t *tmp = matDup(trans);
      matMul(tmp,rep->Gen[i]);
      matMul(tmp,bi);
      mrAddGenerator(result, tmp, 0);
   }
   matFree(bi);
   return result;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
