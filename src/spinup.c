/* ============================= C MeatAxe ==================================
   File:        $Id: spinup.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Spin-up.
   --------------------------------------------------------------------------
   (C) Copyright 1998-2002 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"
#include <string.h>

MTX_DEFINE_FILE_INFO


/**
 ** @defgroup spinup Spin-up and Split
 ** @{
 ** @par Spin-Up
 ** Given a matrix representation and a seed vector v, the spin-up algorithm
 ** calculates the submodule generated by the seed vector, i.e., the smallest 
 ** subspace containing v which is invariant under the generators.
 ** SpinUp() can handle multiple seed vectors, search for cyclic vectors
 ** generating the whole space, and generate seed vectors as linear combinations
 ** of a given basis.
 **
 ** @par Spin-up Scripts
 ** When spinning up a seed vector, you can record the operations performed
 ** by the algorithm in a spin-up script. This script can then be fed into
 ** SpinUpWithScript() to repeat the procedure with a different seed vector 
 ** and different generators.
 **
 ** @par Standard Basis
 ** Normally, the basis vectors computed during the spin-up process are chosen
 ** randomly. However, the spin-up algorithm can be used in "standard basis" mode.
 ** In this mode, the result is invariant under a change of basis.
 ** More precisely, if a given seed vector v and generators g<sub>1</sub>,...g<sub>n</sub>
 ** produce the basis (b<sub>1</sub>,...b_<sub>m</sub>), and A is a nonsingular matrix,
 ** then vA and A<sup>-1</sup>g<sub>1</sub>A,...A<sup>-1</sup>g<sub>n</sub>A produce the basis
 ** (b<sub>1</sub>A,...b<sub>m</sub>A).
 **/


/**
 ** @class SpinUpInfo_t
 ** Spin-up Parameters.
 ** The SpinUpInfo_t data structure is used to pass additional parameters
 ** to the spin-up algorithm, and to return extended results to the caller.
 ** @c Result is set by SpinUp() to report the success of the spin-up.
 ** Possible values are 0 (successful), 1 (not found), and -1 (error).
 **/

/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

static int Dim = 0;			/* Dimension of the whole space */
static int *Piv = NULL;			/* Pivot table */
static const Matrix_t *Seed = NULL;	/* Seed space */
static Matrix_t *Span = NULL;		/* Span */
static int SpanDim = 0;			/* Dimension of the span */
static int Flags = 0;			/* Flags */
static int NGen;			/* Number of generators */
static const Matrix_t **Gen;		/* Generators */
static const Perm_t **GenP;		/* Generators (permutation mode) */
static Matrix_t *StdSpan = NULL;	/* Span */
static long *Script = NULL;

#define OPVEC(i) Script[2*(i)]
#define OPGEN(i) Script[2*(i)+1]


/* --------------------------------------------------------------------------
   Spin1() - Spin up one seed vector

   Arguments:
     <seed>: The seed vector.
     <seedno>: The seed vector number.

   Description:
     This function spins up the next seed vector.

   Return:
     0 = success, 1 = no success, -1 = error
   -------------------------------------------------------------------------- */

static int Spin1(PTR seed, int seedno, const SpinUpInfo_t *info)
{
    int igen;	    /* Generator to apply next */
    PTR get;	    /* Vector to map */
    PTR put;	    /* New vectors go here */
    PTR stdget, stdput;
    FEL f;
    int iget;
    int maxdim;
    int num_tries = 0;	/* Number of multiplications */

    /* Set maxdim: If the user specified a limit, use it. Otherwise use <Dim>.
       ---------------------------------------------------------------------- */
    if (info != NULL && info->MaxSubspaceDimension > 0)
	maxdim = info->MaxSubspaceDimension + 1;
    else
	maxdim = Dim + 1;


    /* If we are not in 'combine' mode, start with an empty space.
       In 'combine' mode, keep what has been found so far.
       ------------------------------------------------------------ */
    if ((Flags & SF_MODE_MASK) != SF_COMBINE)
    {
	SpanDim = 0;
    }

    /* Initialize <get> and <put> to point to the first free row in <Span>
       ------------------------------------------------------------------- */
    get = MatGetPtr(Span,SpanDim);
    iget = SpanDim;
    put = get;
    if (Flags & SF_STD)
    {
	stdget = MatGetPtr(StdSpan,SpanDim);
	stdput = stdget;
    }

    /* Copy the seed vector to <put>, extending the pivot table
       -------------------------------------------------------- */
    FfCopyRow(put,seed);
    FfCleanRow(put,Span->Data,SpanDim,Piv);
    if (Script != NULL)
    {
	OPVEC(SpanDim) = seedno;
	OPGEN(SpanDim) = -1;	/* -1 means it's a seed vector */
    }
    if ((Piv[SpanDim] = FfFindPivot(put,&f)) >= 0)
    {
	++SpanDim;
	FfStepPtr(&put);
	if (Flags & SF_STD)
	{
    	    FfCopyRow(stdput,seed);
	    FfStepPtr(&stdput);
	}
    }

    /* Spin up
       ------- */
    igen = 0;
    while (   get != put && SpanDim < Dim && SpanDim < maxdim && NGen > 0
	   && (   info == NULL || info->MaxTries <= 0
	       || num_tries < info->MaxTries))
    {
	/* Apply the next generator to <get>. IF this was 
	   the last generator, advance <get> to the next vector.
	   ----------------------------------------------------- */
	if (Flags & SF_STD)
	{
	    if (Gen != NULL)
		FfMapRow(stdget,Gen[igen]->Data,Dim,stdput);
	    else
		FfPermRow(stdget,GenP[igen]->Data,stdput);
	    FfCopyRow(put,stdput);
	}
	else
	{
	    if (Gen != NULL)
		FfMapRow(get,Gen[igen]->Data,Dim,put);
	    else
		FfPermRow(get,GenP[igen]->Data,put);
	}
    	if (Script != NULL)
    	{
    	    OPVEC(SpanDim) = iget;
	    OPGEN(SpanDim) = igen;
	}
        if (++igen >= NGen)
	{
	    ++num_tries;
	    igen = 0;
	    FfStepPtr(&get);
	    FfStepPtr(&stdget);
	    ++iget;
	    if (MSG4 || (MSG3 && SpanDim % 50 == 0))
            	printf("SpinUp(): dim=%d, stack=%d\n", SpanDim,SpanDim - iget);
	}

	/* Clean the result with the existing basis
	   ---------------------------------------- */
	FfCleanRow(put,Span->Data,SpanDim,Piv);
	if ((Piv[SpanDim] = FfFindPivot(put,&f)) >= 0)
	{
	    num_tries = 0;
	    ++SpanDim;
	    FfStepPtr(&put);
	    FfStepPtr(&stdput);
	}
    }

    /* Return success code, depending on the mode and the result 
       of the spin up 
       --------------------------------------------------------- */
    MESSAGE(2,("SpinUp(): sub=%d, quot=%d\n",SpanDim,Dim-SpanDim));
    switch (Flags & SF_MODE_MASK)
    {
	case SF_SUB:
	    return (SpanDim > 0 && SpanDim < Dim && SpanDim < maxdim) ?
		    0 : 1;
	case SF_CYCLIC:
	case SF_COMBINE:
	    return SpanDim < Dim ? 1 : 0;
    }
    MTX_ERROR1("Invalid search mode %d",Flags & SF_MODE_MASK);
    return -1;
}



/* --------------------------------------------------------------------------
   CheckArgs0() - Check arguments (common to matrix and permutation mode)
   -------------------------------------------------------------------------- */

static int CheckArgs0(const Matrix_t *seed, int flags)
{
    if (!MatIsValid(seed))
	return -1;
    if (seed->Nor < 1) 
    {
	MTX_ERROR("Empty seed space");
	return -1;
    }
    if (flags == -1)
    {
	MTX_ERROR("Invalid flags");
	return -1;
    }
    return 0;
}




/* --------------------------------------------------------------------------
   CheckArgs() - Check arguments
   -------------------------------------------------------------------------- */

static int CheckArgs(const Matrix_t *seed, const MatRep_t *rep, int flags)
{
    if (CheckArgs0(seed,flags) != 0)
	return -1;
    if (!MrIsValid(rep))
	return -1;
    if (rep->NGen < 0)
    {
	MTX_ERROR1("Invalid number of generators (%d)",rep->NGen);
	return -1;
    }
    if (rep->NGen > 0)
    {
	if (rep->Gen[0]->Noc != seed->Noc || rep->Gen[0]->Field != seed->Field)
	{
	    MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
	    return -1;
	}
    }
    return 0;
}



/* --------------------------------------------------------------------------
   CheckArgsP() - Check arguments (permutation mode)
   -------------------------------------------------------------------------- */

static int CheckArgsP(const Matrix_t *seed, int ngen, const Perm_t **gen, int flags)
{
    int i;
    if (CheckArgs0(seed,flags) != 0)
	return -1;
    if (ngen < 0)
    {
	MTX_ERROR1("Invalid number of generators (%d)",ngen);
	return -1;
    }
    for (i = 0; i < ngen; ++i)
    {
	if (!PermIsValid(gen[i]))
	    return -1;
	if (gen[i]->Degree != seed->Noc)
	{
	    MTX_ERROR3("Gen=%d, seed=%d: %E",gen[i]->Degree,seed->Noc,
		MTX_ERR_INCOMPAT);
	    return -1;
	}
    }
    return 0;
}



/* --------------------------------------------------------------------------
   Init0() - Initialize for spin-up (common for matrix and permutation mode)
   -------------------------------------------------------------------------- */

static int Init0(const Matrix_t *seed, int flags, IntMatrix_t **script,
     const SpinUpInfo_t *info)
{
    int ws_size;
    FfSetField(seed->Field);
    FfSetNoc(seed->Noc);
    Flags = flags;
    Dim = seed->Noc;
    Seed = seed;

    /* Allocate workspace, assuming the worst case.
       -------------------------------------------- */
    if (info && info->MaxSubspaceDimension > 0)
    	ws_size = info->MaxSubspaceDimension + 1;
    else
	ws_size = Dim + 1;
    Span = MatAlloc(seed->Field,ws_size,Dim);
    if (Span == NULL)
    {
	MTX_ERROR("Cannot allocate result buffer");
	return -1;
    }


    SpanDim = 0;
    Piv = NREALLOC(Piv,int,Dim+2);
    if (Piv == NULL)
    {
	MTX_ERROR("Cannot allocate pivot table");
	return -1;
    }
    if (script != NULL)
    {
	if (*script != NULL && 
	    ((*script)->Noc != 2 || (*script)->Nor < Dim + 1))
	{
	    ImatFree(*script);
	    *script = NULL;
	}
	if (*script == NULL)
	    *script = ImatAlloc(Dim+1,2);
	if (*script == NULL)
	{
	    MTX_ERROR("Cannot allocate script");
	    return -1;
	}
	Script = (*script)->Data;
    }
    else
	Script = NULL;	/* No script */

    if (Flags & SF_STD)
    {
	StdSpan = MatAlloc(seed->Field,Dim+1,Dim);
	if (StdSpan == NULL)
	{
	    MatFree(Span);
	    MTX_ERROR("Cannot allocate workspace");
	    return -1;
	}
    }
    return 0;
}





/* --------------------------------------------------------------------------
   Init() - Initialize for spin-up
   -------------------------------------------------------------------------- */

static int Init(const Matrix_t *seed, const MatRep_t *rep, int flags, 
    IntMatrix_t **script, const SpinUpInfo_t *info)
{
    if (Init0(seed,flags,script,info) != 0)
	return -1;
    Gen =  (const Matrix_t **) rep->Gen;
    GenP = NULL;
    NGen = rep->NGen;
    return 0;
}




/* --------------------------------------------------------------------------
   InitP() - Initialize for spin-up (permutation mode)
   -------------------------------------------------------------------------- */

static int InitP(const Matrix_t *seed, int ngen, const Perm_t **gen, int flags, 
    IntMatrix_t **script, SpinUpInfo_t *info)
{
    if (Init0(seed,flags,script,info) != 0)
	return -1;
    GenP = gen;
    Gen = NULL;
    NGen = ngen;
    return 0;
}


/* --------------------------------------------------------------------------
   DoSpinup() - Spin up

   Return:
     0 = Ok, 1 = Not found, -1 = Error
   -------------------------------------------------------------------------- */

int DoSpinup(const SpinUpInfo_t *info)
{
    long n;
    int result;
    PTR seed;

    switch (Flags & SF_SEED_MASK)
    {
    	case SF_FIRST:
	    /* Try the first seed vector only
	       ------------------------------ */
	    return Spin1(Seed->Data,1,info);

    	case SF_EACH:
	    /* Try each seed vector until successful
	       ------------------------------------- */
	    for (seed = Seed->Data, n = 1; n <= Seed->Nor; FfStepPtr(&seed),++n)
	    {
	    	if (Spin1(seed,n,info) == 0)
		    return 0;
	    }
	    return 1;

    	case SF_MAKE:
	    /* Try each 1-dimensional subspace until successfull
	       ------------------------------------------------- */
	    if ((seed = FfAlloc(1)) == NULL) 
	    {
		MTX_ERROR("Cannot allocate seed vector");
		return -1;
	    }
	    result = 1;
	    for (n = 0; result == 1 && (n = MakeSeedVector(Seed,n,seed)) > 0; )
	    {
	    	if (Spin1(seed,n,info) == 0)
		    result = 0;
	    }
	    SysFree(seed);
	    return result;
   
	default:
	    MTX_ERROR("Invalid seed mode");
    }
    return -1;
}





static Matrix_t *DoIt(IntMatrix_t **script, SpinUpInfo_t *info)
{
    int rc = DoSpinup(info);
    if (info != NULL)
	info->Result = rc;

    if (rc < 0)
    {
	MatFree(Span);
	if (Flags & SF_STD)
	    MatFree(StdSpan);
	MTX_ERROR("Spin-up failed");
	return NULL;
    }

    /* Adjust the result size
       ----------------------- */
    if (Flags & SF_STD)
    {
	MatFree(Span);
	Span = StdSpan;
    }
    else
    {
	MatEchelonize(Span);
    }
    Span->Nor = SpanDim;
    Span->Data = (PTR) SysRealloc(Span->Data,FfCurrentRowSize * SpanDim);
    if (script != NULL)
    {
	(*script)->Data = NREALLOC((*script)->Data,long,2 * SpanDim);
	(*script)->Nor = SpanDim;
    }

    return Span;
}






/**
 ** Spin up.
 ** This function calculates the submodule generated by one or more "seed"
 ** vectors under the action of a set of matrices. @a seed must be a matrix
 ** with the same number of columns as the generators and any number of rows.
 ** Of course, all matrices, generators and seed, must be over the same field.
 **
 ** The spinup mode and various options are controlled by two arguments,
 ** @a flags and @a info. @a flags must be a combination of the following
 ** values:
 ** - @c SF_FIRST: Only the first row of @a seed is taken as seed vector.
 ** - @c SF_EACH:  Each row of @a seed is taken as seed vector.
 ** - @c SF_MAKE:  One vector from each 1-dimensional subspace of the row space of @a seed
 **     is taken as seed vector.
 ** - @c SF_SUB: Find a submodule: spin up seed vectors one-by-one until a seed
 **     vector generates a proper submodule.
 ** - @c SF_CYCLIC: Find a cyclic vector: spin up vectors one-by-one until a seed
 **     vector generates the whole space.
 ** - @c SF_COMBINE: Calculate the submodule generated by the set of all seed vectors.
 **     This ist typically used with @c SF_EACH to calculate the submodule
 **     generate by the row space of @a seed.
 ** - @c SF_STD: Create the standard basis. This increases both computation time and
 **     memory usage.
 **
 ** The seed modes, @c SF_FIRST, @c SF_EACH and @c SF_MAKE, and the
 ** search modes, @c SF_SUB, @c SF_CYCLIC, @c SF_COMBINE, are mutually exclusive.
 ** If, in mode @c SF_SUB or @c SF_CYCLIC, no seed vector generates a proper
 ** submodule or the whole space, respectively, this is not considered
 ** an error, and the return value is not @c NULL. The rows of the matrix
 ** returned by SpinUp() always form a basis of an invariant subspace,
 ** but you must examine the number of rows of that matrix to find out if it
 ** is a proper subspace, or null, or the whole space.

 ** The subspace returned by SpinUp() is always in echelon form, if @c SF_STD
 ** is not used. With @c SF_STD however, the subspace is not necessarily in
 ** echelon form.

 ** SpinUp() can record the operations that led to the invariant subspace
 ** in a "spin-up script". You can use the script as input to
 ** SpinUpWithScript() to repeat the spin-up with a different seed
 ** vector. Typically, a spin-up script is created together with @em SF_STD,
 ** and then used to reconstruct the standard basis in a different
 ** representation.
 ** In order to create a spin-up script, @a script must point to a variable
 ** of type IntMatrix_t*. This variable must either be 0 or contain
 ** a valid pointer. In the second case, the buffer pointed to by @a script
 ** is first deallocated before a new script is created. After SpinUp()
 ** returns, the variable contains a pointer to the script. If no spinup
 ** script is needed, pass 0 as 4th parameter.

 ** The format of the spinup script is a matrix with 2 columns and one row
 ** for each basis vector. A row (n,-1) means that the corresponding basis
 ** vector is the n-th seed vector. Seed vector numbers start from 1.
 ** An entry (n,g) with g≥0 means
 ** that the corresponding basis vector was obtained by multiplying the n-th
 ** basis vector by the g-th generator.  Basis vector number and generator
 ** number start from 0.

 ** Additional parameters can be passed via the @a info argument. To be
 ** compatible with future versions of SpinUpInfo_t, you should always
 ** initialize the parameter structure with SpinUpinfoInit().
 **
 ** @param seed Matrix with seed vectors.
 ** @param rep Pointer to a MatRep_t structure with generators.
 ** @param flags Flags, a combination of @c SF_XXXX constants (see description).
 ** @param script Pointer to a variable where the spinup script will be stored (see
    description). May be 0.
 ** @param info Pointer to a data structure with additional parameters, or 0.
 ** @return Span of the seed vector(s) under the action of the generators, or 0
    on error.
 **/

Matrix_t *SpinUp(const Matrix_t *seed, const MatRep_t *rep, int flags,
    IntMatrix_t **script, SpinUpInfo_t *info)
{
    if (CheckArgs(seed,rep,flags) != 0)
    {
	MTX_ERROR1("%E",MTX_ERR_BADARG);
	return NULL;
    }
    if (Init(seed,rep,flags,script,info) != 0)
    {
	MTX_ERROR("Initialization failed");
	return NULL;
    }
    return DoIt(script,info);
}



/**
 ** Spin Up With Permutations.
 ** This function works like Spinup() but expects permutations instead of matrices
 ** for the generators.
 **/

Matrix_t *SpinUpWithPermutations(const Matrix_t *seed, int ngen, 
    const Perm_t **gen, int flags, IntMatrix_t **script, SpinUpInfo_t *info)

{
    if (CheckArgsP(seed, ngen, gen, flags) != 0)
    {
	MTX_ERROR1("%E",MTX_ERR_BADARG);
	return NULL;
    }
    if (InitP(seed,ngen,gen,flags,script,info) != 0)
    {
	MTX_ERROR("Initialization failed");
	return NULL;
    }
    return DoIt(script,info);
}







/**
 ** Initialize spin-up parameters.
 ** @param info Pointer to parameter structure.
 ** @return 0 on success, -1 on error.
 **/

int SpinUpInfoInit(SpinUpInfo_t *info)
{
    memset(info,0,sizeof(*info));
    info->MaxSubspaceDimension = -1;
    return 0;
}



/**
 ** @}
 **/
