////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Calculate homogenous parts of a module
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Makes the standard basis for each basis vector of the peak word kernel.
/// Returns the list of standard bases or NULL on error.

static Matrix_t** MkStdBases(Matrix_t* NPW, MatRep_t* M, const IntMatrix_t* op)
{
   const int num_seed = NPW->nor;
   Matrix_t** V = NALLOC(Matrix_t*, num_seed);
   for (int i = 0; i < num_seed; ++i) {
      Matrix_t* seed = matDupRows(NPW, i, 1);
      V[i] = SpinUpWithScript(seed, M, op);
      matFree(seed);
      if (V[i] == NULL) {
         mtxAbort(MTX_HERE, "SpinUpWithScript() failed for vector %d", i);
      }
   }
   return V;
}
    
/// @addtogroup algo
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Homogeneous part of a module.
/// @param m The module, M.
/// @param s An irreducible constituent of M.
/// @param npw Null-space of the peak word.
/// @param op Spin-up script for the standard basis of S.
/// @param dimends Dimension of the endomorphism ring of S.
/// @return A basis of the S-homogeneous part of $M$, or NULL on error.

Matrix_t* HomogeneousPart(
   MatRep_t* m, MatRep_t* s, Matrix_t* npw,
   const IntMatrix_t* op, int dimends)
{
   Matrix_t
   * mat,
   * a, * b;
   PTR row, vec, basptr;

   const int fl = s->Gen[0]->field;
   const uint32_t Sdim = s->Gen[0]->nor;
   const uint32_t Mdim = m->Gen[0]->nor;
   const uint32_t nulldim = npw->nor;
   MTX_ASSERT(op->nor == Sdim);
   Matrix_t** V = MkStdBases(npw, m, op);

   // Make the system of equations. Nullspace of A gives the desired submodule.
   const uint32_t len = Mdim * m->NGen * Sdim;  // The number of equations
   MTX_LOG2("HomogeneousPart(): len=%lu", (unsigned long)len);
   Matrix_t* A = matAlloc(fl, nulldim, len);
   for (int i = 0; i < m->NGen; i++) {
      uint32_t colin = i * Mdim * Sdim;     // the first place in a row that is to fill
      MTX_LOG2("colin=%d, nulldim=%d, Sdim=%d", colin, nulldim, Sdim);
      for (uint32_t j = 0; j < nulldim; j++) {
         PTR matptr = matGetPtr(A, j);
         int u;
         a = matDup(V[j]);
         b = matDup(s->Gen[i]);
         matMul(a, m->Gen[i]);                  /* the equations that describe  */
         matMul(b, V[j]);                       /* that a vector in the null-   */
         matMulScalar(b, ffNeg(FF_ONE));        /* space is the first element   */
         matAdd(a, b);                          /* of a standard basis of a     */
         /* module isomorphic to S       */
         uint32_t col = colin;                /* the index of the temporary column  of A */
         for (u = 0; u < Sdim; ++u) {
            PTR v = matGetPtr(a, u);
            int t;
            for (t = 0; t < Mdim; col++, t++) {
               FEL f = ffExtract(v, t);
               ffInsert(matptr, col, f);
            }
         }
         matFree(a);
         matFree(b);
      }
   }

   MTX_LOG2("Equation system is %lux%lu", (unsigned long)A->nor, (unsigned long)A->noc);
   Matrix_t* gensys = matNullSpace__(A);     /* module-generating system for the S-part */

   // spin up the basis of the whole S-part of M
   MTX_ASSERT(Sdim % dimends == 0);
   const uint32_t dim = gensys->nor * (Sdim / dimends);
   MTX_ASSERT(dim % Sdim == 0);
   uint32_t nr = dim / Sdim;
   Matrix_t* bas = matAlloc(fl, dim, Mdim);
   basptr = bas->data;
   vec = gensys->data;
   for (uint32_t i = 1; i <= gensys->nor; i++, ffStepPtr(&vec, nulldim)) {
      Matrix_t* seed = matAlloc(fl, 1, Mdim);
      for (uint32_t j = 0; j < nulldim; j++) {
         FEL f = ffExtract(vec, j);
         mat = matDup(V[j]);
         row = mat->data;
         ffMulRow(row, f, Mdim);
         ffAddRow(seed->data, row, Mdim);
         matFree(mat);
      }
      Matrix_t* base = matDup(bas); // basis for the whole S-isomorphic part of M
      matEchelonize(base);

      if (!IsSubspace(seed, base, 0)) {
         PTR v;
         // Copy into bas the standard basis for one S-isomorphic submodule of M.
         nr--;
         Matrix_t* sum = matAlloc(fl, Sdim, Mdim);
         for (uint32_t j = 0; j < nulldim; j++) {        // calculates the SB to u
            FEL f;
            MTX_ASSERT(j < gensys->noc);
            f = ffExtract(vec, j);
            matAddMul(sum, V[j], f);
         }
         v = sum->data;
         for (uint32_t j = 0; j < Sdim; j++) {
            ffCopyRow(basptr, v, Mdim);
            ffStepPtr(&basptr, Mdim);
            ffStepPtr(&v, Mdim);
         }
         matFree(sum);
      }
      matFree(base);
      matFree(seed);
      // if (!nr) break;	// the whole basis is found
   }

   matFree(gensys);
   for (uint32_t i = 0; i < npw->nor; ++i) {
      matFree(V[i]);
   }
   sysFree(V);

   return bas;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
