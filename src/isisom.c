////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Isomorphism test
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"




////////////////////////////////////////////////////////////////////////////////////////////////////
/// @addtosection algo
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

static void checkArgs(int ngen, Matrix_t  **gen1, const CfInfo *info1, Matrix_t **gen2, int use_pw)
{
   int j;

   /* Check arguments
      --------------- */
   MTX_ASSERT(ngen > 0);
   for (j = 0; j < ngen; ++j)
   {
      matValidate(MTX_HERE, gen1[j]);
      matValidate(MTX_HERE, gen2[j]);
      if (gen1[j]->nor != gen1[j]->noc)
      {
         mtxAbort(MTX_HERE,"gen1[%d]: Matrix not square",j);
      }
      if (gen2[j]->nor != gen2[j]->noc)
      {
         mtxAbort(MTX_HERE,"gen2[%d]: Matrix not square",j);
      }
      if (gen1[j]->field != gen1[0]->field || gen1[j]->nor != gen1[0]->nor)
      {
         mtxAbort(MTX_HERE,"gen1[%d]: Incompatible matrix",j);
      }
      if (gen2[j]->field != gen1[0]->field)
      {
         mtxAbort(MTX_HERE,"gen2[%d]: Incompatible matrix",j);
      }
   }

   if (info1->dim != gen1[0]->nor)
   {
      mtxAbort(MTX_HERE,"Inconsistent cfinfo data");
   }
   if (use_pw && info1->peakWord == 0)
   {
      mtxAbort(MTX_HERE,"No peak word available");
   }
   if (!use_pw && info1->idWord == 0)
   {
      mtxAbort(MTX_HERE,"No id word available");
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Compare representations.
/// This function decides if two irreducible representations are isomorphic.
/// |rep1| and |rep2| must be two matrix representations over the same field,
/// and with the same number of generators. Furthermore,
/// to compare the representations, the function needs an identifying word
/// for the first representation, i.e., the fields |info1->idWord|, 
/// |info1->idPol| and |info1->spl| must be set, and the generators in 
/// |rep1| must be in standard basis with respect to the identifying word.
/// If |use_pw| is nonzero, the peak word is used instead of the idword.
/// In this case, |rep1| must of course be in standard basis with respect 
/// to the peak word.
///
/// If the representations are isomorphic, and |trans| is not |NULL|, the
/// basis transformation which makes the second representation identical
/// to the first is stored into |*trans|. To be more precise, if $g_i$ is 
/// the representation of the i-th generator in representation |rep1|,
/// h<sub>i</sub> in representation @a rep2, and T the 
/// matrix returned in @a trans, then $T h<sub>i</sub> T<sup>-1</sup>=g<sub>i</sub>.
/// @param rep1 The first representation.
/// @param info1 Pointer to an CfInfo_t structure for the first representation.
/// @param rep2 The second representation.
/// @param trans Buffer for basis transformation matrix, or NULL.
/// @param use_pw If different from zero, use peak word instead of the identifying word.
/// @return 1 if the representations are isomorphic, 0 otherwise

int IsIsomorphic(const MatRep_t *rep1, const CfInfo *info1,
    		 const MatRep_t *rep2, Matrix_t  **trans, int use_pw)
{
    int j;
    WgData_t *wg;
    Matrix_t  *word, *m, *seed, *b, *bi;
    int result;

    checkArgs(rep1->NGen,rep1->Gen,info1,rep2->Gen,use_pw);
 
    /* Check if the dimensions are equal
       --------------------------------- */
    if (rep1->Gen[0]->nor != rep2->Gen[0]->nor)
	return 0;

    /* Make the idword on representation 2
       ----------------------------------- */
    wg = wgAlloc(rep2);
    word = wgMakeWord(wg,use_pw ? info1->peakWord : info1->idWord);
    m = matInsert(word,use_pw ? info1->peakPol : info1->idPol);
    matFree(word);
    wgFree(wg);
    seed = matNullSpace__(m);
    if (seed->nor != info1->spl)
    {	
	matFree(seed);
	return 0;
    }

    /* Make the standard basis
       ----------------------- */
    b = SpinUp(seed,rep2,SF_FIRST|SF_CYCLIC|SF_STD,NULL,NULL);
    matFree(seed);
    if (b->nor != b->noc)
    {
	matFree(b);
	return 0;
    }
    bi = matInverse(b);

    /* Compare generators
       ------------------ */
    for (j = 0, result = 0; result == 0 && j < rep2->NGen; ++j)
    {
	Matrix_t *g = matDup(b);
	matMul(g,rep2->Gen[j]);
	matMul(g,bi);
	if (matCompare(g,rep1->Gen[j]) != 0)
	    result = 1;
	matFree(g);
    }

    /* Clean up 
       -------- */
    if (trans != NULL && result == 0)
	*trans = b;
    else
	matFree(b);
    matFree(bi);

    return (result == 0);
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
