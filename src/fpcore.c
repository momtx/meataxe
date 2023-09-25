////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Basic factored polynomial functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


#define FP_MAGIC 0x17B69244

/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @class FPoly_t
/// A Factored Polynomial.
/// This structure contains a polynomial which is split into factors. The factors
/// need not be irreducible.

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check a factored polynomial.
/// @param p The polynomial.
/// @return 1 if @em p is a valid factored polynomial, 0 otherwise.

int fpIsValid(const FPoly_t *p)
{
   int i;
   if (p == NULL) {
      return 0;
   }
   if ((p->typeId != FP_MAGIC) || (p->NFactors < 0) || (p->BufSize < p->NFactors)) {
      return 0;
   }
   if ((p->Factor == NULL) || (p->Mult == NULL)) {
      return 0;
   }
   for (i = 0; i < p->NFactors; ++i) {
      if (!polIsValid(p->Factor[i]))
         return 0;
      if (p->Mult[i] < 0) {
         return 0;
      }
      if ((i > 0) && (p->Factor[i]->field != p->Factor[0]->field)) {
         return 0;
      }
   }
   return 1;
}

void fpValidate(const struct MtxSourceLocation* src, const FPoly_t *p)
{
   int i;
   if (p == NULL) {
      mtxAbort(src,"NULL polynomial");
   }
   if ((p->typeId != FP_MAGIC) || (p->NFactors < 0) || (p->BufSize < p->NFactors)) {
      mtxAbort(src,"Invalid FPoly_t: Magic=%d, NFactors=%d, MaxLen=%d",
                 (int)p->typeId,p->NFactors,p->BufSize);
   }
   if ((p->Factor == NULL) || (p->Mult == NULL)) {
      mtxAbort(src,"Invalid FPoly_t: Factor:%s, Mult:%s",
                 p->Factor == 0 ? "NULL" : "ok",
                 p->Mult == 0 ? "NULL" : "ok");
   }
   for (i = 0; i < p->NFactors; ++i) {
      polValidate(src, p->Factor[i]);
      if (p->Mult[i] < 0) {
         mtxAbort(src,"Invalid multiplicity %d",p->Mult[i]);
      }
      if ((i > 0) && (p->Factor[i]->field != p->Factor[0]->field)) {
         mtxAbort(src,"Factors over different fields");
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Allocate a factored polynomial.
/// This function creates a new Fpoly_t structure.
/// The new polynomial is empty, i.e., it has no factors.
/// @return Pointer to the new FPoly_t structure.

FPoly_t *fpAlloc()
{
   FPoly_t* x = ALLOC(FPoly_t);
   x->BufSize = 5;
   x->Factor = NALLOC(Poly_t *,x->BufSize);
   x->Mult = NALLOC(int,x->BufSize);
   x->NFactors = 0;
   x->typeId = FP_MAGIC;
   return x;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Free a factored polynomial.
/// @return 0 on success, -1 on error.
/// @see FPoly_t FpAlloc

int fpFree(FPoly_t *x)
{
   int i;

   /* Check the argument
      ------------------ */
   fpValidate(MTX_HERE, x);

   /* Free all factors
      ---------------- */
   for (i = 0; i < x->NFactors; ++i) {
      polFree(x->Factor[i]);
   }

   /* Free the <FPoly_t> structure
      ---------------------------- */
   sysFree(x->Factor);
   sysFree(x->Mult);
   memset(x,0,sizeof(FPoly_t));
   sysFree(x);
   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
