////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Polynomial g.c.d.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


////////////////////////////////////////////////////////////////////////////////////////////////////

static void normalize(Poly_t *p, FEL f)
{
   register FEL *buf;
   register int i;

   if (f == FF_ONE) {
      return;
   }
   for (buf = p->data, i = (int) p->degree; i >= 0; --i, ++buf) {
      if (*buf != FF_ZERO) {
         *buf = ffDiv(*buf,f);
      }
   }
}


/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Greatest common divisor of two polynomials.
/// This function calculates the gratest common divisor of two poynomials.
/// The polynomials must be over the same field, and at least one of
/// them must be different from zero.
/// Unlike most polynomial functions, polGcd() normalizes the result,
/// i.e., the leading coefficient of the g.c.d., is always one.
/// @see polGcdEx()
/// @param a First polynomial.
/// @param b Second polynomial.
/// @return Greatest common divisor of @em a and @em b, 0 on error.

Poly_t *polGcd(const Poly_t *a, const Poly_t *b)
{
   Poly_t *p,*q,*tmp;
   FEL f;

   // check arguments
   polValidate(MTX_HERE, a);
   polValidate(MTX_HERE, b);
   if (a->field != b->field) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }

   // Handle special cases
   if (a->degree == -1) {
      if (b->degree == -1) {
         mtxAbort(MTX_HERE,"%s",MTX_ERR_DIV0);
         return NULL;
      }
      return polDup(b);
   } else if (b->degree == -1) {
      return polDup(a);
   }

   // calculate g.c.d.
   ffSetField(a->field);
   if (a->degree < b->degree) {
      p = polDup(b);
      q = polDup(a);
   } else {
      p = polDup(a);
      q = polDup(b);
   }
   while (q->degree >= 0) {
      if (polMod(p,q) == NULL) {
         return NULL;
      }
      tmp = p;
      p = q;
      q = tmp;
   }
   polFree(q);

   // normalize the result
   if ((f = p->data[p->degree]) != FF_ONE) {
      normalize(p,f);
   }

   return p;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Greatest common divisor of two polynomials.
/// Given two polynomials a and b, this function calculates the
/// greatest common divisor g=gcd(a,b) and two polynomials p, q
/// such that g=pa+qb. Both @em a and @em b must be nonzero. The leading
/// coefficient of g is always one.
///
/// @em result must be a pointer to an array of three <tt>Poly_t *</tt> elements.
/// If the function is successful, pointers to g, p, and q have been stored in
/// result[0], result[1], and result[2], respectively. The caller is
/// responsible for destroying these polynomials when they are no longer
/// needed. In particular, you must not use the same |result| buffer again
/// with %polGcdEx() before the contents have been freed.
/// @see polGcd()
/// @param a First polynomial.
/// @param b Second polynomial.
/// @param result Result buffer for the g.c.d. and coefficients (see below).
/// @return 0 on sucess, -1 on error.

int polGcdEx(const Poly_t *a, const Poly_t *b, Poly_t **result)
{
   Poly_t *p, *q, *pa, *pb, *qa, *qb;
   FEL f;
   int alessb;

   // check arguments
   polValidate(MTX_HERE, a);
   polValidate(MTX_HERE, b);
   if (result == NULL) {
      return mtxAbort(MTX_HERE,"result: %s",MTX_ERR_BADARG), -1;
   }
   if (a->field != b->field) {
      return mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT), -1;
   }

   alessb = a->degree < b->degree;
   p = polDup(alessb ? b : a);
   q = polDup(alessb ? a : b);
   pb = polAlloc(a->field,alessb ? 0 : -1);
   qa = polAlloc(a->field,alessb ? 0 : -1);
   pa = polAlloc(a->field,alessb ? -1 : 0);
   qb = polAlloc(a->field,alessb ? -1 : 0);

   while (q->degree >= 0) {
      int i;
      Poly_t *tmp, *atmp, *btmp;
      Poly_t *m = polDivMod(p,q);
      tmp = p;
      p = q;
      q = tmp;

      atmp = polDup(qa);
      btmp = polDup(qb);
      for (i = 0; i <= m->degree; ++i) {
         m->data[i] = ffSub(FF_ZERO,m->data[i]);
      }
      polMul(atmp,m);
      polMul(btmp,m);
      polAdd(atmp,pa);
      polAdd(btmp,pb);

      polFree(pa);
      polFree(pb);

      pa = qa;
      pb = qb;

      qa = atmp;
      qb = btmp;
      polFree(m);
   }

   // normalize the result
   if ((f = p->data[p->degree]) != FF_ONE) {
      normalize(p,f);
      normalize(pa,f);
      normalize(pb,f);
   }

   result[0] = p;
   result[1] = pa;
   result[2] = pb;
   polFree(q);
   polFree(qa);
   polFree(qb);

   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
