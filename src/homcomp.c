////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Calculate homogenous parts of a module
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"




////////////////////////////////////////////////////////////////////////////////////////////////////

/// Makes the standard basis for each basis vector of the peak word kernel.
/// Returns the list of standard bases or NULL on error.

static Matrix_t **MkStdBasis(Matrix_t *NPW, MatRep_t *M, const IntMatrix_t *op)
{
   const int num_seed = NPW->Nor;

   Matrix_t **V = NALLOC(Matrix_t *,num_seed);
   int ok = (V != NULL);

   for (int i = 0; ok && i < num_seed; ++i)
   {
    	Matrix_t *seed = matCutRows(NPW,i,1);
    	V[i] = SpinUpWithScript(seed,M,op);
	matFree(seed);
	if (V[i] == NULL) {
	    mtxAbort(MTX_HERE,"SpinUpWithScript() failed for vector %d",i);
            return V;
        }
   }

   // Clean up
   if (!ok && V != NULL) {
      for (int i = 0; i < num_seed; ++i) {
         if (V[i])
            matFree(V[i]);
      }
      sysFree(V);
      V = NULL;
   }

   return V;
}
    






    
/// @addtosection algo
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Homogeneous part of a module.
/// @param m The module, M.
/// @param s An irreducible constituent of M.
/// @param npw Null-space of the peak word.
/// @param op Spin-up script for the standard basis of S.
/// @param dimends Dimension of the endomorphism ring of S.
/// @return A basis of the S-homogeneous part of $M$, or NULL on error.

Matrix_t *HomogeneousPart(MatRep_t *m, MatRep_t *s, Matrix_t *npw, 
                          const IntMatrix_t *op, int dimends)
{
   Matrix_t
	    **V,	    /*The standard basis for the given basis NPW */
	    *A,		    /* The matrix, the nullspace of which gives the
			       desired submodule */
	    *sum,	    /* The standard basis for one submodule of M
			       that is isomorphic to S */
	    *bas, *base,    /* basis for the whole S-isomorphic part of M */
	    *gensys,	    /* a module-generating system for the S-part  */
	    *mat, *seed, 
	    *a, *b;
   int i, nr, dim, Sdim, nulldim, Mdim, fl, col, len;
   PTR row, vec, basptr;

    fl = s->Gen[0]->Field;
    Sdim = s->Gen[0]->Nor;
    Mdim = m->Gen[0]->Nor;
    nulldim = npw->Nor;
    MTX_ASSERT(op->Nor == Sdim, NULL);
    V = MkStdBasis(npw,m,op);

    /* Make the system of equations.
       ----------------------------- */
    len = Mdim * m->NGen * Sdim;		/* The number of equations */
    MESSAGE(3,("HomogeneousPart(): len=%d\n",len));
    if ((A = matAlloc(fl, nulldim, len)) == NULL)
    {
	mtxAbort(MTX_HERE,"Cannot allocate buffer");
	return NULL;
    }
    for (i = 0; i < m->NGen; i++)
    {
        int colin, j;
	colin = i * Mdim * Sdim;   /* the first place in a row that is to fill */
	MESSAGE(3,("colin=%d, nulldim=%d, Sdim=%d\n",colin,nulldim,Sdim));
	for (j=0; j<nulldim; j++)
	{
	    PTR matptr = matGetPtr(A,j);
	    int u;
	    a = matDup(V[j]);
	    b = matDup(s->Gen[i]);
	    matMul(a,m->Gen[i]);		/* the equations that describe  */
	    matMul(b,V[j]);			/* that a vector in the null-   */
	    matMulScalar(b,ffNeg(FF_ONE));	/* space is the first element   */
	    matAdd(a, b);			/* of a standard basis of a     */ 
					/* module isomorphic to S       */
	    col = colin;     	     /* the index of the temporary column  of A */
	    for (u=0; u < Sdim; ++u)
	    {
		PTR v = matGetPtr(a,u);
		int t;
		for (t = 0; t < Mdim; col++, t++)
		{
		    FEL f = ffExtract(v,t);
		    ffInsert(matptr,col,f);
		}
	    }
	    matFree(a);
	    matFree(b);
	}
    }

    MESSAGE(2,("Equation system is %dx%d\n",A->Nor,A->Noc));
    gensys = matNullSpace__(A);    /* module-generating system for the S-part */

/* spins up the basis of the whole S-part of M
   ------------------------------------------- */
    MTX_ASSERT(Sdim % dimends == 0, NULL);
    dim = gensys->Nor * (Sdim/dimends);
    MTX_ASSERT(dim % Sdim == 0, NULL);
    nr = dim/Sdim;		
    if ((bas = matAlloc(fl, dim, Mdim)) == NULL)
    {
	mtxAbort(MTX_HERE,"Cannot allocate buffer");
	return NULL;
    }
    basptr = bas->Data;
    vec = gensys->Data;
    for (i = 1; i <= gensys->Nor; i++, ffStepPtr(&vec, nulldim))
    {
    	int j;
	if ((seed = matAlloc (fl,1,Mdim)) == NULL)
	{
	    mtxAbort(MTX_HERE,"Cannot allocate buffer");
	    return NULL;
	}
	for (j = 0; j < nulldim; j++)
	{
	    FEL f;
	    f = ffExtract(vec, j);
	    mat = matDup(V[j]);
	    row = mat->Data;
	    ffMulRow(row, f, Mdim);
	    ffAddRow(seed->Data, row, Mdim);
	    matFree(mat);
	}
	base = matDup(bas);
	matEchelonize(base);

	if (!IsSubspace(seed,base,0))
	{
	    PTR v;
	    /* Copy into bas the standard basis for 
	       one S-isomorphic submodule of M.
               -------------------------------- */
	    matFree(seed);
	    nr--;
	    if((sum = matAlloc(fl, Sdim, Mdim)) == NULL)
	    {
		mtxAbort(MTX_HERE,"Cannot allocate buffer");
		return NULL;
	    }
	    for (j = 0; j < nulldim; j++)	/* calculates the SB to u */
	    {
		FEL f;
		MTX_ASSERT(j < gensys->Noc, NULL);
		f = ffExtract(vec,j);
		matAddMul(sum,V[j],f);
	    }
	    v = sum->Data;
	    for (j = 0; j < Sdim; j++)
	    {
		ffCopyRow(basptr, v, Mdim);
		ffStepPtr(&basptr, Mdim);
		ffStepPtr(&v, Mdim);
	    }
	    matFree(sum);
	}
	matFree(base); 
/* if (!nr) break;	*/	/* the whole basis is found */
    }

    return bas;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
