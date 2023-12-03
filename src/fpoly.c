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
/// A factored polynomial.
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
   if ((p->typeId != FP_MAGIC) || (p->nFactors < 0) || (p->bufSize < p->nFactors)) {
      return 0;
   }
   if ((p->factor == NULL) || (p->mult == NULL)) {
      return 0;
   }
   for (i = 0; i < p->nFactors; ++i) {
      if (!polIsValid(p->factor[i]))
         return 0;
      if (p->factor[i]->field != p->field) {
         return 0;
      }
      if (p->mult[i] < 0) {
         return 0;
      }
      if ((i > 0) && (p->factor[i]->field != p->factor[0]->field)) {
         return 0;
      }
   }
   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void fpValidate(const struct MtxSourceLocation* src, const FPoly_t* p)
{
   int i;
   if (p == NULL) {
      mtxAbort(src, "NULL polynomial");
   }
   if ((p->typeId != FP_MAGIC) || (p->nFactors < 0) || (p->bufSize < p->nFactors)) {
      mtxAbort(src, "Invalid FPoly_t: Magic=%d, nFactors=%d, MaxLen=%d",
         (int)p->typeId, p->nFactors, p->bufSize);
   }
   if ((p->factor == NULL) || (p->mult == NULL)) {
      mtxAbort(src, "Invalid FPoly_t: factor:%s, mult:%s",
         p->factor == 0 ? "NULL" : "ok",
         p->mult == 0 ? "NULL" : "ok");
   }
   for (i = 0; i < p->nFactors; ++i) {
      polValidate(src, p->factor[i]);
      if (p->factor[i]->field != p->field) {
         mtxAbort(src, "Invalid FPoly_t: Inconsistent field orders (%lu vs %lu)",
            (unsigned long) p->field,
            (unsigned long) p->factor[i]->field);
      }

      if (p->mult[i] < 0) {
         mtxAbort(src, "Invalid FPoly_t: Invalid multiplicity %d", p->mult[i]);
      }
      if ((i > 0) && (p->factor[i]->field != p->factor[0]->field)) {
         mtxAbort(src, "Invalid FPoly_t: Factors over different fields");
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Allocate a factored polynomial.
/// This function creates a new Fpoly_t structure.
/// The new polynomial is empty, i.e., it has no factors.
/// @return Pointer to the new FPoly_t structure.

FPoly_t *fpAlloc(uint32_t field)
{
   FPoly_t* x = ALLOC(FPoly_t);
   x->field = field;
   x->bufSize = 5;
   x->factor = NALLOC(Poly_t *,x->bufSize);
   x->mult = NALLOC(int,x->bufSize);
   x->nFactors = 0;
   x->typeId = FP_MAGIC;
   return x;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Creates a copy of a factored polynomial.

FPoly_t *fpDup(const FPoly_t *src)
{
   fpValidate(MTX_HERE, src);

   // Copy the factors
   Poly_t **new_factor = NALLOC(Poly_t *,src->nFactors);
   int* new_mult = NALLOC(int,src->nFactors);
   for (int i = 0; i < src->nFactors; ++i) {
      new_mult[i] = src->mult[i];
      new_factor[i] = polDup(src->factor[i]);
   }

   // Create a new factored polynomial
   FPoly_t *x = fpAlloc(src->field);
   sysFree(x->factor);
   sysFree(x->mult);
   x->factor = new_factor;
   x->mult = new_mult;
   x->bufSize = x->nFactors = src->nFactors;
   return x;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Compare two factored polynomials.

int fpCompare(const FPoly_t* a, const FPoly_t* b)
{
   if (a->field > b->field) {
      return 1;
   }
   if (a->field < b->field) {
      return -1;
   }
   int i;
   for (i = 0; i < a->nFactors && i < b->nFactors; ++i)
   {
      int cmp = polCompare(a->factor[i], b->factor[i]);
      if (cmp != 0) {
         return cmp;
      }
      if (a->mult[i] > b->mult[i]) {
         return 1;
      }
      if (a->mult[i] < b->mult[i]) {
         return -1;
      }
   }
   if (i < a->nFactors) {
      return 1;
   }
   if (i < b->nFactors) {
      return -1;
   }
   return 0;
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
   for (i = 0; i < x->nFactors; ++i) {
      polFree(x->factor[i]);
   }

   /* Free the <FPoly_t> structure
      ---------------------------- */
   sysFree(x->factor);
   sysFree(x->mult);
   memset(x,0,sizeof(FPoly_t));
   sysFree(x);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply with an irreducible polynomial.
/// This function multiplies a factored polynomial with the power of an
/// an irreducible factor. It is not checked that @em src is irreducible.
/// @see fpMul()
/// @param dest Factored polynomial to modify.
/// @param src Irreducible polynomial.
/// @param pwr Power of the irreducible polynomial.
/// @return The function returns @em dest or 0 on error.

FPoly_t* fpMulP(FPoly_t* dest, const Poly_t* src, int pwr)
{
   int i;
   int cmp = 0;

   // Check the arguments
   polValidate(MTX_HERE, src);
   fpValidate(MTX_HERE, dest);
   if (src->field != dest->field) {
      mtxAbort(MTX_HERE, "Inconsistent fields (%lu vs %lu)",
         (unsigned long) src->field, (unsigned long) dest->field);
      return NULL;
   }
   if (pwr <= 0) {
      mtxAbort(MTX_HERE, "pwr=%d: %s", pwr, MTX_ERR_BADARG);
      return NULL;
   }

   // Find the insert position
   for (i = 0;
        i < dest->nFactors && (cmp = polCompare(dest->factor[i], src)) < 0;
        ++i) {}

   // Extend the buffer, if necessary
   if ((i >= dest->nFactors) || (cmp != 0)) {
      int k;
      if (dest->nFactors >= dest->bufSize) {
         int newsize = dest->bufSize + 5;
         Poly_t** x = NREALLOC(dest->factor, Poly_t*, newsize);
         int* e = NREALLOC(dest->mult, int, newsize);
         dest->factor = x;
         dest->mult = e;
         dest->bufSize = newsize;
      }

      // Make room for the new factor
      for (k = dest->nFactors; k > i; --k) {
         dest->factor[k] = dest->factor[k - 1];
         dest->mult[k] = dest->mult[k - 1];
      }
      ++dest->nFactors;

      // Insert new factor
      dest->factor[i] = polDup(src);
      dest->mult[i] = pwr;
      if (dest->factor[i] == NULL) {
         mtxAbort(MTX_HERE, "Cannot copy polynomial");
         return NULL;
      }
   }
   else {
      dest->mult[i] += pwr;
   }
   return dest;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Multiply factored polynomials.
/// Multiplies @em dest by @em src. The previous content of @em dest is lost.
/// @see fpMulP()
/// @param dest Factored polynomial to modify.
/// @param src Factored polynomial.
/// @return The function returns |dest| or |NULL| on error.

FPoly_t *fpMul(FPoly_t *dest, const FPoly_t *src)
{
   int i;

   fpValidate(MTX_HERE, src);
   fpValidate(MTX_HERE, dest);

   for (i = 0; i < src->nFactors; ++i) {
      fpMulP(dest,src->factor[i],src->mult[i]);
   }
   return dest;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Format a factored polynomial

void fpFormat(StrBuffer* sb, const FPoly_t *p)
{
   fpValidate(MTX_HERE, p);
   if (p->nFactors == 0) {
      sbAppend(sb, "1");
      return;
   }

   for (int i = 0; i < p->nFactors; ++i) {
      int e = p->mult[i];
      if (i > 0) { sbAppend(sb, " * "); }
      sbAppend(sb, "(");
      polFormat(sb, p->factor[i]);
      if (e > 1) {
         sbPrintf(sb, ")^%d",e);
      } else {
         sbAppend(sb,")");
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////


/// Print a factored polynomial.
/// This function prints a factored polynomial to the standard output.
/// If @em name is not 0, "name=" is printed before the polynomial.
/// @param name Name of the polynomial or 0.
/// @param p Pointer to the factored polynomial.

void fpPrint(const char *name, const FPoly_t *p)
{
   fpValidate(MTX_HERE, p);
   StrBuffer* sb = sbAlloc(100);
   if (name != NULL) {
      sbPrintf(sb, "%s =",name);
   }
   fpFormat(sb, p);
   fputs(sbData(sb), stdout);
   sbFree(sb);
   if (name != NULL) {
      printf("\n");
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

char* fpToEphemeralString(const FPoly_t *p)
{
   StrBuffer* sb = sbAlloc(100);
   fpFormat(sb, p);
   return sbToEphemeralString(sb);
}


/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
