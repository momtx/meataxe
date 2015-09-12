////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Polynomial g.c.d.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

////////////////////////////////////////////////////////////////////////////////////////////////////

static void normalize(Poly_t *p, FEL f)
{
   register FEL *buf;
   register int i;

   if (f == FF_ONE) {
      return;
   }
   for (buf = p->Data, i = (int) p->Degree; i >= 0; --i, ++buf) {
      if (*buf != FF_ZERO) {
         *buf = FfDiv(*buf,f);
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
/// Unlike most polynomial functions, PolGcd() normalizes the result,
/// i.e., the leading coefficient of the g.c.d., is always one.
/// @see PolGcdEx()
/// @param a First polynomial.
/// @param b Second polynomial.
/// @return Greatest common divisor of @em a and @em b, 0 on error.

Poly_t *PolGcd(const Poly_t *a, const Poly_t *b)
{
   Poly_t *p,*q,*tmp;
   FEL f;

   // check arguments
   if (!PolIsValid(a) || !PolIsValid(b)) {
      return NULL;
   }
   if (a->Field != b->Field) {
      MTX_ERROR1("%E",MTX_ERR_INCOMPAT);
      return NULL;
   }

   // Handle special cases
   if (a->Degree == -1) {
      if (b->Degree == -1) {
         MTX_ERROR1("%E",MTX_ERR_DIV0);
         return NULL;
      }
      return PolDup(b);
   } else if (b->Degree == -1) {
      return PolDup(a);
   }

   // calculate g.c.d.
   FfSetField(a->Field);
   if (a->Degree < b->Degree) {
      p = PolDup(b);
      q = PolDup(a);
   } else {
      p = PolDup(a);
      q = PolDup(b);
   }
   while (q->Degree >= 0) {
      if (PolMod(p,q) == NULL) {
         return NULL;
      }
      tmp = p;
      p = q;
      q = tmp;
   }
   PolFree(q);

   // normalize the result
   if ((f = p->Data[p->Degree]) != FF_ONE) {
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
/// with %PolGcdEx() before the contents have been freed.
/// @see PolGcd()
/// @param a First polynomial.
/// @param b Second polynomial.
/// @param result Result buffer for the g.c.d. and coefficients (see below).
/// @return 0 on sucess, -1 on error.

int PolGcdEx(const Poly_t *a, const Poly_t *b, Poly_t **result)
{
   Poly_t *p, *q, *pa, *pb, *qa, *qb;
   FEL f;
   int alessb;

   // check arguments
   if (!PolIsValid(a) || !PolIsValid(b)) {
      return -1;
   }
   if (result == NULL) {
      return MTX_ERROR1("result: %E",MTX_ERR_BADARG), -1;
   }
   if (a->Field != b->Field) {
      return MTX_ERROR1("%E",MTX_ERR_INCOMPAT), -1;
   }

   alessb = a->Degree < b->Degree;
   p = PolDup(alessb ? b : a);
   q = PolDup(alessb ? a : b);
   pb = PolAlloc(a->Field,alessb ? 0 : -1);
   qa = PolAlloc(a->Field,alessb ? 0 : -1);
   pa = PolAlloc(a->Field,alessb ? -1 : 0);
   qb = PolAlloc(a->Field,alessb ? -1 : 0);

   while (q->Degree >= 0) {
      int i;
      Poly_t *tmp, *atmp, *btmp;
      Poly_t *m = PolDivMod(p,q);
      tmp = p;
      p = q;
      q = tmp;

      atmp = PolDup(qa);
      btmp = PolDup(qb);
      for (i = 0; i <= m->Degree; ++i) {
         m->Data[i] = FfSub(FF_ZERO,m->Data[i]);
      }
      PolMul(atmp,m);
      PolMul(btmp,m);
      PolAdd(atmp,pa);
      PolAdd(btmp,pb);

      PolFree(pa);
      PolFree(pb);

      pa = qa;
      pb = qb;

      qa = atmp;
      qb = btmp;
      PolFree(m);
   }

   // normalize the result
   if ((f = p->Data[p->Degree]) != FF_ONE) {
      normalize(p,f);
      normalize(pa,f);
      normalize(pb,f);
   }

   result[0] = p;
   result[1] = pa;
   result[2] = pb;
   PolFree(q);
   PolFree(qa);
   PolFree(qb);

   return 0;
}


/// @}
