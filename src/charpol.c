////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Characteristic and minimal polynomial of a matrix
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

/// @defgroup charpol Characteristic and Minimal Polynomials
/// @{

static const uint32_t CHARPOL_MAGIC = 0xb92bbdc5;

////////////////////////////////////////////////////////////////////////////////////////////////////

void charpolValidate(const struct MtxSourceLocation* src, Charpol_t* cp)
{
   if (cp == NULL || cp->magic != CHARPOL_MAGIC)
      mtxAbort(src,"Invalid charpol state");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Starts computation of the characteristic or minimal polynomial of a matrix.
/// This function returns a computation state object which can be passed to @ref charpolFactor
/// to calculate the characteristic or minimal polynomial of a matrix.
///
/// @param matrix The matrix.
/// @param mode Decides which polynomial (minimal or characteristic) to compute.

Charpol_t* charpolAlloc(const Matrix_t* matrix, enum CharpolMode mode, long seed)
{
   matValidate(MTX_HERE, matrix);
   if (matrix->Nor != matrix->Noc) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTSQUARE);
   }
   Charpol_t* state = ALLOC(Charpol_t);
   
   state->magic = CHARPOL_MAGIC;
   state->mode = mode;
   state->fl = matrix->Field;
   state->vsDim = matrix->Nor;
   ffSetField(state->fl);
   state->mat = ffAlloc(state->vsDim, state->vsDim);
   state->A = ffAlloc(state->vsDim + 1, state->vsDim);
   state->B = ffAlloc(state->vsDim + 1, state->vsDim);
   state->piv = NALLOC(long,state->vsDim + 2);
   state->ispiv = NALLOC(char,state->vsDim + 2);
       
   // TODO: provide option to suppress copy and use the original matrix.
   // This should always done in charPol() and may be useful for charPolFactor().
   memcpy(state->mat,matrix->Data,ffSize(state->vsDim,state->vsDim));
   memset(state->ispiv,0,(size_t)(state->vsDim + 2));

   state->seed = seed < state->vsDim ? seed : 0;

   if (mode == PM_MINPOL)
      state->partialMinPol = polAlloc(state->fl, 0);

   return state;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Destroys a polynomial computation state.
/// Call this function when the polynomial computation is finished or cancelled. All memory
/// associated with the computation is released, and the @a state pointer becomes invalid.

void charpolFree(struct CharpolState* state)
{
   charpolValidate(MTX_HERE, state);
   sysFree(state->mat);
   sysFree(state->A);
   sysFree(state->B);
   sysFree(state->piv);
   sysFree(state->ispiv);
   if (state->partialMinPol != NULL) polFree(state->partialMinPol);
   memset(state, 0, sizeof(*state));
   sysFree(state);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Make polynomial for the latest cyclic subspace

static Poly_t *mkpoly(struct CharpolState* state)
{
   Poly_t *pol = polAlloc(state->fl,state->n);
   PTR x = ffGetPtr(state->B,state->n, state->vsDim);
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
   uint32_t pv;
   FEL f;

   PTR a = ffGetPtr(state->A, state->dim, state->vsDim);
   PTR b = state->B;
   ffMulRow(b,FF_ZERO, state->vsDim);
   state->n = 0;
   while ((pv = ffFindPivot(a,&f, state->vsDim)) != MTX_NVAL) {
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
      ffStepPtr(&a, state->vsDim);
      ffMapRow(x,state->mat,state->vsDim,state->vsDim,a);
      y = b;
      ffStepPtr(&b, state->vsDim);
      ffMulRow(b,FF_ZERO, state->vsDim);
      for (uint32_t k = 0; k < state->vsDim - 1; ++k) {
         ffInsert(b,k + 1,ffExtract(y,k));
      }

      /* Clean with existing basis vectors
         --------------------------------- */
      x = state->A;
      y = state->B;
      for (uint32_t k = 0; k < state->dim + state->n; ++k) {
         f = ffDiv(ffExtract(a,state->piv[k]),ffExtract(x,state->piv[k]));
         ffAddMulRow(a,x,ffNeg(f), state->vsDim);
         if (k >= state->dim) {
            ffAddMulRow(b,y,ffNeg(f), state->vsDim);
            ffStepPtr(&y, state->vsDim);
         }
         ffStepPtr(&x, state->vsDim);
      }
   }
   state->dim += state->n;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static Poly_t *charpolFactor_(struct CharpolState* state)
{
   MTX_ASSERT(state != NULL);
   MTX_ASSERT(state->mat);
   if (state->dim >= state->vsDim)      // nothing left to do
   {
      return NULL;
   }

   // Prepare the next seed vector
   ffSetField(state->fl);
   /*    seed = ffGetPtr(A,state->dim,state->nor);*/
   PTR seed = (PTR)((char *)state->A + ffSize(state->dim, state->vsDim));
   int i;
   if (state->dim == 0) {
      i = state->seed;
   } else {
      for (i = 0; i < state->vsDim && state->ispiv[i] != 0; ++i);
   }
   ffMulRow(seed,FF_ZERO, state->vsDim);
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

/// Computes one factor of the characteristic or minimal polynomial of a matrix.
///
/// The function needs a computation state for the matrix which must have been created with
/// @ref charpolAlloc. Each call of @c charpolFactor for the same state returns a new factor of the
/// characteristic or minimal polynomial. If the polynomial is complete, the function returns NULL.
///
/// The factors returned by @c charpolFactor are in general reducible. If you need the
/// characteristic polynomial in fully factored form, use @ref charpol.

Poly_t *charpolFactor(Charpol_t* state)
{
   charpolValidate(MTX_HERE, state);
   return charpolFactor_(state);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the charateristic polynomial of a matrix in fully factorized form.

FPoly_t *charpol(const Matrix_t *mat)
{
   matValidate(MTX_HERE, mat);
   Charpol_t* state = charpolAlloc(mat, PM_CHARPOL, 0);

   FPoly_t *cpol = fpAlloc();
   Poly_t *p;
   for (p = charpolFactor_(state); p != NULL; p = charpolFactor_(state)) {
      FPoly_t *factors = Factorization(p);
      polFree(p);
      fpMul(cpol,factors);
      fpFree(factors);
   }
   charpolFree(state);
   return cpol;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Computes one factor of the minimal polynomial of a matrix.
///
/// The function needs a computation state for the matrix which must have been created with
/// @ref charpolAlloc using @ref PM_MINPOL mode. Each call of @c charpolFactor for the same state
/// returns a new factor of the minimal polynomial. If the polynomial is complete, the function
/// returns NULL
///
/// The factors returned by @c minpolFactor are in general reducible. If you need the
/// characteristic polynomial in fully factored form, use @ref charpol.

Poly_t *minpolFactor(Charpol_t* state)
{
   charpolValidate(MTX_HERE, state);
   MTX_ASSERT(state->mode == PM_MINPOL);
   return charpolFactor_(state);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the minimal polynomial of a matrix in fully factorized form.

FPoly_t *minpol(const Matrix_t *mat)
{
   matValidate(MTX_HERE, mat);
   Charpol_t* state = charpolAlloc(mat, PM_MINPOL, 0);
   FPoly_t *mp = fpAlloc();
   while (1)
   {
      Poly_t* f = charpolFactor_(state);
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
