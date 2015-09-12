////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Calculate homogenous parts of a module
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>


MTX_DEFINE_FILE_INFO


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Makes the standard basis for each basis vector of the peak word kernel.

static Matrix_t **MkStdBasis(Matrix_t *NPW, MatRep_t *M, const IntMatrix_t *op)
{
    Matrix_t **V;
    int num_seed = NPW->Nor;
    int i;

    if ((V = NALLOC(Matrix_t *,num_seed)) == NULL)
    {
	MTX_ERROR("Cannot allocate buffer");
	return NULL;
    }
    for (i = 0; i < num_seed; ++i)
    {
    	Matrix_t *seed = MatCutRows(NPW,i,1);
    	V[i] = SpinUpWithScript(seed,M,op);
	if (V[i] == NULL)
	    MTX_ERROR1("SpinUpWithScript() failed for vector %d",i);
	MatFree(seed);
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
    MTX_VERIFY(op->Nor == Sdim);
    V = MkStdBasis(npw,m,op);

    /* Make the system of equations.
       ----------------------------- */
    len = Mdim * m->NGen * Sdim;		/* The number of equations */
    MESSAGE(3,("HomogeneousPart(): len=%d\n",len));
    if ((A = MatAlloc(fl, nulldim, len)) == NULL)
    {
	MTX_ERROR("Cannot allocate buffer");
	return NULL;
    }
    for (i = 0; i < m->NGen; i++)
    {
        int colin, j;
	colin = i * Mdim * Sdim;   /* the first place in a row that is to fill */
	MESSAGE(3,("colin=%d, nulldim=%d, Sdim=%d\n",colin,nulldim,Sdim));
	for (j=0; j<nulldim; j++)
	{
	    PTR matptr = MatGetPtr(A,j);
	    int u;
	    a = MatDup(V[j]);
	    b = MatDup(s->Gen[i]);
	    MatMul(a,m->Gen[i]);		/* the equations that describe  */
	    MatMul(b,V[j]);			/* that a vector in the null-   */
	    MatMulScalar(b,FfNeg(FF_ONE));	/* space is the first element   */
	    MatAdd(a, b);			/* of a standard basis of a     */ 
					/* module isomorphic to S       */
	    col = colin;     	     /* the index of the temporary column  of A */
	    FfSetNoc(len);
	    for (u=0; u < Sdim; ++u)
	    {
		PTR v = MatGetPtr(a,u);
		int t;
		for (t = 0; t < Mdim; col++, t++)
		{
		    FEL f = FfExtract(v,t);
		    FfInsert(matptr,col,f);
		}
	    }
	    MatFree(a);
	    MatFree(b);
	    FfSetNoc(len);   
	}
    }

    MESSAGE(2,("Equation system is %dx%d\n",A->Nor,A->Noc));
    gensys = MatNullSpace__(A);    /* module-generating system for the S-part */

/* spins up the basis of the whole S-part of M
   ------------------------------------------- */
    MTX_VERIFY(Sdim % dimends == 0);
    dim = gensys->Nor * (Sdim/dimends);
    MTX_VERIFY(dim % Sdim == 0);
    nr = dim/Sdim;		
    if ((bas = MatAlloc(fl, dim, Mdim)) == NULL)
    {
	MTX_ERROR("Cannot allocate buffer");
	return NULL;
    }
    basptr = bas->Data;
    vec = gensys->Data;
    FfSetNoc(nulldim);
    for (i = 1; i <= gensys->Nor; i++, FfStepPtr(&vec))
    {
    	int j;
	if ((seed = MatAlloc (fl,1,Mdim)) == NULL)
	{
	    MTX_ERROR("Cannot allocate buffer");
	    return NULL;
	}
	for (j = 0; j < nulldim; j++)
	{
	    FEL f;
	    f = FfExtract(vec, j);
	    mat = MatDup(V[j]);
	    row = mat->Data;
	    FfSetNoc(Mdim);
	    FfMulRow(row, f);
	    FfAddRow(seed->Data, row);
	    MatFree(mat);
	}
	base = MatDup(bas);
	MatEchelonize(base);

	if (!IsSubspace(seed,base,0))
	{
	    PTR v;
	    /* Copy into bas the standard basis for 
	       one S-isomorphic submodule of M.
               -------------------------------- */
	    MatFree(seed);
	    nr--;
	    if((sum = MatAlloc(fl, Sdim, Mdim)) == NULL)
	    {
		MTX_ERROR("Cannot allocate buffer");
		return NULL;
	    }
	    for (j = 0; j < nulldim; j++)	/* calculates the SB to u */
	    {
		FEL f;
		MTX_VERIFY(j < gensys->Noc);
		f = FfExtract(vec,j);
		MatAddMul(sum,V[j],f);
	    }
	    v = sum->Data;
	    FfSetNoc(Mdim);
	    for (j = 0; j < Sdim; j++)
	    {
		FfCopyRow(basptr, v);
		FfStepPtr(&basptr);
		FfStepPtr(&v);
	    }
	    MatFree(sum);
	}
	MatFree(base); 
/* if (!nr) break;	*/	/* the whole basis is found */
	FfSetNoc(nulldim);
    }

    return bas;
}


/// @}
