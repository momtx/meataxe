////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Characteristic and minimal polynomial of a matrix
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

enum PolyMode { PM_CHARPOL, PM_MINPOL };

/// @defgroup charpol Characteristic and Minimal Polynomials
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Seed for Characteristic Polynomial.
/// This variable is used by charPolFactor() to select the first
/// seed vector. By default, CharPolSeed has the value 0, i.e., the
/// first seed vector is (1,0,...,0). Assigning the value 1 selects
/// the start vector (0,1,...,0) in all subsequent calls to
/// charPolFactor().
/// If CharPolSeed is out of bounds, charPolFactor() will reset it to 0.

long CharPolSeed = 0;           /* Seed */

/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

struct CharpolState 
{
   enum PolyMode mode;

   long fl;      // field order
   long nor;     // vector space dimension  TODO: rename to vsDim
   long *piv;    // Pivot table
   char *ispiv;  // Pivot flags
   PTR mat;      // The matrix
   PTR A;        // Work space (for spin-up)
   PTR B;        // Work space II (coefficients)
   long dim;     // Dimension reached so far (sum of cyclic subspace dimensions) TODO: rename to currentDim
   long n;       // Dimension of cyclic subspace
   long seed;    // Number if seed vector for the first cyclic subspace

   // Minimal polynomial for the current subspace dimension. NULL in characteristic polynomial mode.
   Poly_t* partialMinPol;  
};

// will be removed for MT support
static struct CharpolState state__;

// will be rewritten for MT support
static void destroyState(struct CharpolState* state)
{
   if (state->mat == NULL)
      return;

   sysFree(state->mat);
   sysFree(state->A);
   sysFree(state->B);
   sysFree(state->piv);
   sysFree(state->ispiv);
   if (state->partialMinPol != NULL) polFree(state->partialMinPol);

   memset(state, 0, sizeof(*state));
}

// will be rewritten for MT support
static struct CharpolState* createState(const Matrix_t* matrix, enum PolyMode mode)
{
   matValidate(MTX_HERE, matrix);
   if (matrix->Nor != matrix->Noc) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTSQUARE);
   }
   
   if (state__.mat != NULL) destroyState(&state__);
   struct CharpolState* state = &state__;
   memset(state, 0, sizeof(*state));

   state->mode = mode;
   state->fl = matrix->Field;
   state->nor = matrix->Nor;
   ffSetField(state->fl);
   state->mat = ffAlloc(state->nor, state->nor);
   state->A = ffAlloc(state->nor + 1, state->nor);
   state->B = ffAlloc(state->nor + 1, state->nor);
   state->piv = NALLOC(long,state->nor + 2);
   state->ispiv = NALLOC(char,state->nor + 2);
       
   // TODO: provide option to suppress copy and use the original matrix.
   // This should always done in charPol() and may be useful for charPolFactor().
   memcpy(state->mat,matrix->Data,ffSize(state->nor,state->nor));
   memset(state->ispiv,0,(size_t)(state->nor + 2));

   state->seed = (CharPolSeed >= 0 && CharPolSeed < state->nor) ? CharPolSeed : 0;
   CharPolSeed = state->seed;  // reproduce MeatAxe v2.3 behaviour

   if (mode == PM_MINPOL)
      state->partialMinPol = polAlloc(state->fl, 0);

   return state;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static struct CharpolState* getState(enum PolyMode mode)
{
   MTX_ASSERT(state__.mat != NULL, NULL);
   MTX_ASSERT(state__.mode == mode, NULL);
   return &state__;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Make polynomial for the latest cyclic subspace

static Poly_t *mkpoly(struct CharpolState* state)
{
   Poly_t *pol = polAlloc(state->fl,state->n);
   PTR x = ffGetPtr(state->B,state->n, state->nor);
   for (int k = 0; k < state->n; ++k) {
      pol->Data[k] = ffExtract(x,k);
   }
   pol->Data[state->n] = FF_ONE;
   return pol;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Spin up one cyclic subspace

static void spinup_cyclic(struct CharpolState* state)
{
   long pv, k;
   FEL f;

   PTR a = ffGetPtr(state->A, state->dim, state->nor);
   PTR b = state->B;
   ffMulRow(b,FF_ZERO, state->nor);
   state->n = 0;
   while ((pv = ffFindPivot(a,&f, state->nor)) >= 0) {
      PTR x, y;

      /* Add new vector to basis
         ----------------------- */
      state->piv[state->dim + state->n] = pv;
      state->ispiv[pv] = 1;
      ffInsert(b,state->n,FF_ONE);
      ++state->n;

      /* Calculate the next vector
         ------------------------- */
      x = a;
      ffStepPtr(&a, state->nor);
      ffMapRow(x,state->mat,state->nor,state->nor,a);
      y = b;
      ffStepPtr(&b, state->nor);
      ffMulRow(b,FF_ZERO, state->nor);
      for (k = 0; k < state->nor - 1; ++k) {
         ffInsert(b,k + 1,ffExtract(y,k));
      }

      /* Clean with existing basis vectors
         --------------------------------- */
      x = state->A;
      y = state->B;
      for (k = 0; k < state->dim + state->n; ++k) {
         f = ffDiv(ffExtract(a,state->piv[k]),ffExtract(x,state->piv[k]));
         ffAddMulRow(a,x,ffNeg(f), state->nor);
         if (k >= state->dim) {
            ffAddMulRow(b,y,ffNeg(f), state->nor);
            ffStepPtr(&y, state->nor);
         }
         ffStepPtr(&x, state->nor);
      }
   }
   state->dim += state->n;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the next factor or NULL if the characteristic polynomial is complete.

static Poly_t *CharPolFactor_(struct CharpolState* state)
{
   MTX_ASSERT(state != NULL, NULL);
   MTX_ASSERT(state->mat, NULL);
   if (state->dim >= state->nor)      // nothing left to do
   {
      return NULL;
   }

   // Prepare the next seed vector
   ffSetField(state->fl);
   /*    seed = ffGetPtr(A,state->dim,state->nor);*/
   PTR seed = (PTR)((char *)state->A + ffSize(state->dim, state->nor));
   int i;
   if (state->dim == 0) {
      i = state->seed;
   } else {
      for (i = 0; i < state->nor && state->ispiv[i] != 0; ++i);
   }
   ffMulRow(seed,FF_ZERO, state->nor);
   ffInsert(seed,i,FF_ONE);

   // Spin up and return the polynomial
   spinup_cyclic(state);
   Poly_t* factor = mkpoly(state);

   if (state->mode == PM_MINPOL) {
      Poly_t* gcd = polGcd(factor, state->partialMinPol);
      Poly_t* minpolFactor = polDivMod(factor, gcd);
      polFree(factor);
      polFree(gcd);
      polMul(state->partialMinPol, minpolFactor);
      return minpolFactor;
   }

   return factor;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Characteristic Polynomial.
/// This function returns one factor of the characteristic polynomial of
/// a given matrix. Further calls with a 0 argument return
/// more factors or 0, if there are no more factors.
/// Note that the factors obtained in this way are in general not irreducible.
///
/// Here is how %charPolFactor() works: If @a mat is different from 0,
/// %charPolFactor() initializes its internal data and starts
/// computing one cyclic subspace. The choice of starting vector for this
/// first subspace depends on the global variable CharPolSeed.
/// Usually, this variable has a value of 0, corresponding to the vector
/// (1,0,...,0). Then, the polynomial of the matrix restricted to
/// that cyclic subspace is constructed and returned to the caller.
///
/// If @a mat is 0 on the next call, %charPolFactor() resumes at
/// the point where it returned the last time, calculates the next cyclic
/// subspace and so on, until the complete space is exhausted.
///
/// @attention Since the function uses static variables to store
/// information across multiple calls, your program must not use
/// %charPolFactor() on more than one matrix at the same time.
/// @param mat Pointer to the matrix.
/// @return A factor of the characteristic polynomial, or NULL if there are no more factors.

Poly_t *charPolFactor(const Matrix_t *mat)
{
   struct CharpolState* state =
      (mat != NULL) ? createState(mat, PM_CHARPOL) : getState(PM_CHARPOL);
   return CharPolFactor_(state);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Characteristic Polynomial.
/// This function calculates the characteristic polynomial of a matrix in fully factorized form.
/// The return value is a pointer to a FPoly_t structure containing the irreducible factors of
/// the characteristic polynomial.
///
/// @param mat Pointer to the matrix.
/// @return The characteristic polynomial of @a mat

FPoly_t *charPol(const Matrix_t *mat)
{
   struct CharpolState* state = createState(mat, PM_CHARPOL);

   FPoly_t *cpol = fpAlloc();
   Poly_t *p;

   for (p = CharPolFactor_(state); p != NULL; p = CharPolFactor_(state)) {
      FPoly_t *factors = Factorization(p);
      polFree(p);
      fpMul(cpol,factors);
      fpFree(factors);
   }
   return cpol;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Poly_t* minPolFactor(const Matrix_t *mat)
{
   struct CharpolState* state = mat != NULL ? createState(mat, PM_MINPOL) : getState(PM_MINPOL);
   return CharPolFactor_(state);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

FPoly_t *minPol(const Matrix_t *mat)
{
    matValidate(MTX_HERE, mat);

    FPoly_t *mp = fpAlloc();
    while (1)
    {
        Poly_t* f = minPolFactor(mat);
        if (f == NULL) break;
    	FPoly_t *ff = Factorization(f);
    	polFree(f);
    	fpMul(mp,ff);
    	fpFree(ff);
    }
    return mp;
}
/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
