////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Polynomials over finit fields.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

/// @defgroup poly Polynomials
/// @details
/// The MeatAxe can work with polynomials over a finite field. A polynomial is represented
/// by a Poly_t structure. Each polynomial carries the field order, i.e., you can work
/// with polynomials over different fields on one program. However, this feature is
/// currently of little use since all standard operations only work on polynomials over the
/// same field, and there is no easy way to identify polynomials over a field and its subfields.
///
/// There is a second representation of polynomials as product of factors, see FPoly_t.

/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// @class Poly_t
/// @details
/// Internally, a polynomial of degree n is represented as an array of n+1 field
/// elements (@a data), where <tt>data[i]</tt> is the coefficient of x<sup>i</sup>.
/// The leading coefficient is always non-zero.
/// The zero polynomial has a degree of -1.

////////////////////////////////////////////////////////////////////////////////////////////////////

static void grow(Poly_t *p, int newdeg)
{
   if (p->degree < newdeg) {
      if (p->bufSize < newdeg + 1) {
         p->bufSize = newdeg + 1;
         p->data = NREALLOC(p->data, FEL, p->bufSize);
      }
      for (int32_t i = p->degree + 1; i <= newdeg; ++i) {
         p->data[i] = FF_ZERO;
      }
      p->degree = newdeg;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Add polynomials.
///
/// This function adds @em src to @em dest and returns dest
/// The polynomials must be over the same field.

Poly_t *polAdd(Poly_t *dest, const Poly_t *src)
{
   polValidate(MTX_HERE, src);
   polValidate(MTX_HERE, dest);
   if (dest->field != src->field)
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
   if (src->degree == -1) {
      return dest;      // src = 0 
   }
   ffSetField(src->field);
   grow(dest,src->degree);
   FEL* s = src->data;
   FEL* d = dest->data;
   for (int32_t i = src->degree; i >= 0; --i) {
      *d = ffAdd(*d,*s++);
      ++d;
   }
   polNormalize(dest);
   return dest;
}

// Gives identical results for ZZZ=0 and ZZZ=1.
static inline int ffCompare(FEL a, FEL b)
{
#if MTX_ZZZ == 0
   return a > b ? 1 : ((a == b) ?  0 : -1);
#elif MTX_ZZZ == 1
   int d = ffToInt(a) - ffToInt(b);
   return d > 0 ? 1 : (d == 0 ?  0 : -1);
#else
   #error
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Compare Two poynomials.
/// This function compares two polynomials and returns 0 if the polynomials are equal,
/// -1 if a<b, or 1 if a>b.  The ordering of polynomials is defined as follows.
/// - If a and b are over different fields, the polynomials over
///   the larger field is greater.
/// - Otherwise, if they have different degrees, the polynomial
///   with the higher degree is greater.
/// - If both field and degree are equal, the result of the comparison is 0
///   if the polynomials are equal. Otherwise it is unspecified if the return value
///   is +1 or -1.
/// @param a First polynomial.
/// @param b Second polynomial.
/// @return 1 if a>b, -1 if a<b, 0 if a=b, or -2 on error.

int polCompare(const Poly_t *a, const Poly_t *b)
{
   int i;

   // check arguments
   polValidate(MTX_HERE, a);
   polValidate(MTX_HERE, b);

   // compare fields
   if (a->field > b->field)
      return 1;
   if (a->field < b->field)
      return -1;

   // compare degrees
   if (a->degree > b->degree) {
      return 1;
   }
   if (a->degree < b->degree) {
      return -1;
   }

   // compare coefficients
   for (i = a->degree; i >= 0; --i) {
       int cmp = ffCompare(a->data[i], b->data[i]);
       if (cmp != 0)
	   return cmp;
   }

   // the polynomials are equal
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check a polynomial.
/// This function checks if the argument is a pointer to a valid polynomial. If the polynomial
/// the function returns 1. Otherwise, an error is signalled and, if the error handler does not
/// terminate the program, the function returns 0.
/// @param p The polynomial to check.
/// @return 1 if @em p points to a valid polynomial, 0 otherwise.

int polIsValid(const Poly_t *p)
{
   if (p == NULL || p->typeId != MTX_TYPE_POLYNOMIAL) {
      return 0;
   }
   const int deg = p->degree;
   if (deg < -1 || p->field < 2 || p->data == NULL || p->bufSize < 0) {
      return 0;
   }
   if (deg >= 0 && p->data[deg] == FF_ZERO) {
      return 0;
   }
   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Checks if the given polynomial is valid and aborts the program if the test fails.

void polValidate(const struct MtxSourceLocation* src, const Poly_t* pol)
{
   if (pol == NULL) {
      mtxAbort(src, "NULL polynomial");
   }
   if (pol->typeId != MTX_TYPE_POLYNOMIAL || pol->degree < -1 || pol->field < 2) {
      mtxAbort(src, "Invalid polynomial (typeId=0x%lx, field=%lu, deg=%ld)",
               (unsigned long) pol->typeId, (unsigned long) pol->field, (long) pol->degree);
   }
   if (pol->data == NULL || pol->bufSize < 0) {
      mtxAbort(src, "Invalid polynomial (data=0x%lx, size=%d)",
               (unsigned long)pol->data, pol->bufSize);
   }
   if (pol->degree >= 0 && pol->data[pol->degree] == FF_ZERO) {
      mtxAbort(src, "Invalid polynomial (leading coefficient is zero)");
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Allocate a polynomial
/// This function creates the polynomial p(x)=x^n over the current field.
/// If n is negative, a zero polynomial (degree -1) is created. The return value is a pointer to a
/// newly allocated Poly_t structure. The caller is responsible for releasing memory by calling
/// @ref polFree when the polynomial is no longer needed.
/// @param field Field order.
/// @param degree Degree of the polynomial.
/// @return Pointer to a new Poly_t structure

Poly_t* polAlloc(uint32_t field, int32_t degree)
{
   if (degree < 0) {degree = -1;}
   MTX_ASSERT(degree < (int32_t) 2147483647L);

   ffSetField(field);
   Poly_t* x = (Poly_t*) mmAlloc(MTX_TYPE_POLYNOMIAL, sizeof(Poly_t));
   x->field = field;
   x->degree = degree;
   x->bufSize = degree + 1;
   x->data = NALLOC(FEL, x->bufSize);
   for (int32_t i = 0; i < degree; ++i) {
      x->data[i] = FF_ZERO;
   }
   if (degree >= 0) {
      x->data[degree] = FF_ONE;
   }
   return x;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Free a polynomial
/// This function frees a polynomial data structure and cleans up all internal data.
/// @param x Pointer to the polynomial.
/// @return $0$ on success, $-1$ on error.

void polFree(Poly_t *x)
{
   polValidate(MTX_HERE, x);
   sysFree(x->data);
   x->data = NULL;
   mmFree(x, MTX_TYPE_POLYNOMIAL);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Remove leading zero coefficients.
/// This function makes sure that the leading coefficient of a polynomial is non-zero.
/// It does <em>not</em> necessarily normaliz the leading coefficient to one!

void polNormalize(Poly_t *p)
{
   int i = p->degree;
   while (i >= 0 && p->data[i] == FF_ZERO) {
      --i;
   }
   p->degree = i;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Derive a polynomial.
/// This function derives a polynomial. Note that the derived polynomial is
/// stored in @em pol, replacing the original polynomial. The following piece of
/// code shows how to keep the original polynomial intact while calculating
/// the derivative:
/// @code
/// Poly_t *pol, *der;
/// ...
/// der = polDerive(polDup(pol));
/// @endcode
/// @param pol Pointer to the polynomial.
/// @return @em pol.

Poly_t *polDerive(Poly_t *pol)
{
   int i, maxdeg = -1;
   register FEL *buf;
   FEL f = FF_ZERO;

   // check argument
   polValidate(MTX_HERE, pol);

   buf = pol->data;
   ffSetField(pol->field);
   for (i = 0; i < pol->degree; ++i) {
      f = ffAdd(f,FF_ONE);
      buf[i] = ffMul(buf[i + 1],f);
      if (buf[i] != FF_ZERO) {
         maxdeg = i;
      }
   }
   pol->degree = maxdeg;
   return pol;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Polynomial division.
///
/// This function performs a polynomial division. Given two polynomials a and
/// b≠0 over the same field, %polDivMod() finds two  polynomials q and r such
/// that a=qb+r, and deg(r)<deg(b). See also @ref polMod.
///
/// The quotient q is returned as the function result. This is a newly
/// allocated polynomial. The caller is responsible for deleting the quotient
/// when it no longer needed.
///
/// The remainder r, is stored in @em a and replaces the original value. If you
/// need to preserve the value of @em a you must make a copy using polDup() before
/// calling %polDivMod(). @em b is not changed.
/// @param a First polynomial (numerator) on call, remainder on return.
/// @param b Second polynomial (denominator).
/// @return The quotient or 0 on error.

Poly_t *polDivMod(Poly_t *a, const Poly_t *b)
{
   Poly_t *q;

   // check arguments
   polValidate(MTX_HERE, a);
   polValidate(MTX_HERE, b);
   if (a->field != b->field) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }
   ffSetField(a->field);
   if (b->degree <= -1) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_DIV0);
      return NULL;
   }
   if (a->degree < b->degree) {
      q = polAlloc(a->field,-1);        // trivial case: Quotient = 0
   } else {
      FEL lead = b->data[b->degree];
      int i, k;

      if (lead == FF_ZERO) {
         mtxAbort(MTX_HERE,"%s",MTX_ERR_DIV0);
         return NULL;
      }
      q = polAlloc(ffOrder,a->degree - b->degree);
      for (i = a->degree; i >= b->degree; --i) {
         FEL qq = ffNeg(ffDiv(a->data[i],lead));
         for (k = 0; k <= b->degree; ++k) {
            a->data[i - k] = ffAdd(a->data[i - k],
                                   ffMul(qq,b->data[b->degree - k]));
         }
         q->data[i - b->degree] = ffNeg(qq);
      }
      polNormalize(a);
   }
   return q;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Polynomial division.
/// This function replaces @em a with the remainder of the division of @em a by @em b.
/// @see polDivMod()
/// @param a First polynomial (numerator) on call, remainder on return.
/// @param b Second polynomial (denominator).
/// @return @em a or NULL on error.

Poly_t *polMod(Poly_t *a, const Poly_t *b)
{
   // check arguments
   polValidate(MTX_HERE,a);
   polValidate(MTX_HERE,b);
   if (a->field != b->field) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }

   ffSetField(a->field);
   if (b->degree <= -1) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_DIV0);
      return NULL;
   }
   if (a->degree >= b->degree) {
      FEL lead = b->data[b->degree];
      int i, k;

      if (lead == FF_ZERO) {
         mtxAbort(MTX_HERE,"%s",MTX_ERR_DIV0);
         return NULL;
      }
      for (i = a->degree; i >= b->degree; --i) {
         FEL qq = ffNeg(ffDiv(a->data[i],lead));
         for (k = 0; k <= b->degree; ++k) {
            a->data[i - k] = ffAdd(a->data[i - k], ffMul(qq,b->data[b->degree - k]));
         }
	 MTX_ASSERT(a->data[i] == FF_ZERO);
      }
      polNormalize(a);
   }
   return a;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns a copy of an existing polynomial.

Poly_t *polDup(const Poly_t *p)
{
   polValidate(MTX_HERE,p);
   Poly_t *y = polAlloc(p->field,p->degree);
   memcpy(y->data,p->data,(p->degree + 1) * sizeof(FEL));
   return y;
}

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

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the greatest common divisor of two polynomials.
///
/// The polynomials must be over the same field, and at least one of them must be different from
/// zero. Unlike most polynomial functions, polGcd() normalizes the result,
/// i.e., the leading coefficient of the g.c.d., is always one.
/// @see polGcdEx()

Poly_t *polGcd(const Poly_t *a, const Poly_t *b)
{
   Poly_t *p,*q,*tmp;
   FEL f;

   // check arguments
   polValidate(MTX_HERE, a);
   polValidate(MTX_HERE, b);
   if (a->field != b->field)
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);

   // Handle special cases
   if (a->degree == -1) {
      if (b->degree == -1)
         mtxAbort(MTX_HERE,"%s",MTX_ERR_DIV0);
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

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Multiply polynomials.
///
/// This function multiplies @em dest by @em src and returns @em dest.
/// The polynomials must be over the same field.
/// @param dest Pointer to the first polynomial.
/// @param src Pointer to the second polynomial.
/// @return @em dest, or 0 on error.

Poly_t *polMul(Poly_t *dest, const Poly_t *src)
{
   FEL *x, *y, *d, *s;
   int di, si;
   size_t xdeg;

   // check arguments
   polValidate(MTX_HERE, src);
   polValidate(MTX_HERE, dest);
   if (dest->field != src->field) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
      return NULL;
   }

   // handle special cases: dest = 0, src = 0
   if (dest->degree == -1) {
      return dest;
   }
   if (src->degree == -1) {
      dest->degree = -1;
      return dest;
   }

   d = dest->data;
   s = src->data;
   xdeg = src->degree + dest->degree;
   ffSetField(src->field);

   // allocate result buffer
   x = NALLOC(FEL,xdeg + 1);
   if (x == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate result");
      return NULL;
   }
   for (di = xdeg, y = x; di >= 0; --di) {
      *y++ = FF_ZERO;
   }

   // multiply
   for (di = 0; di <= dest->degree; ++di) {
      for (si = 0; si <= src->degree; ++si) {
         x[si + di] = ffAdd(x[si + di],ffMul(s[si],d[di]));
      }
   }

   // overwrite <dest> with the result
   sysFree(dest->data);
   dest->data = x;
   dest->degree = xdeg;
   dest->bufSize = xdeg + 1;
   return dest;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Format a polynomial
///
/// This function converts a polynomial to a human-readable text form and appends the text to the
/// given string buffer.

void polFormat(StrBuffer_t* sb, const Poly_t *p)
{
   polValidate(MTX_HERE, p);
   ffSetField(p->field);
   if (p->degree == -1) {
      sbPrintf(sb, "0x^0");
   }
   const char* plus = "";
   for (int i = p->degree; i >= 0; i--) {
      if (p->data[i] != FF_ZERO) {
         sbAppend(sb, plus);
         if ((p->data[i] != FF_ONE) || (i == 0)) {
            sbPrintf(sb, "%d",ffToInt(p->data[i]));
         }
         if (i > 1) {
            sbPrintf(sb, "x^%d",i);
         } else if (i == 1) {
            sbAppend(sb, "x");
         }
         plus = "+";
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Print a polynomial
///
/// This function prints a polynomial on the standard output in a human-readable form (see
/// @ref polFormat).
/// If @p name is not 0, the name followed by an equal sign is printed before the polynomial.
/// For example, the statement <tt>polPrint("P",P)</tt> could
/// produce the following output:
/// <pre>
/// P=3x^2+x+1</pre>
///
/// @param name Name to print before the polynomial or 0.
/// @param p Pointer to the polynomial.

void polPrint(char *name, const Poly_t *p)
{
   polValidate(MTX_HERE, p);
   if (name != NULL) {
      printf("%s=",name);
   }
   StrBuffer_t* sb = sbAlloc(30);
   polFormat(sb, p);
   fputs(sbData(sb), stdout);
   sbFree(sb);
   if (name != NULL) {
      printf("\n");
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the human-readable form (as defined by @ref polFormat) of a polynomial in an ephemeral
/// string (see @ref sbToEphemeralString).

char* polToEphemeralString(const Poly_t *p)
{
   StrBuffer_t* sb = sbAlloc(100);
   polFormat(sb, p);
   return sbToEphemeralString(sb);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes a polynomial to a file. See also @ref polSave.

void polWrite(const Poly_t* p, MtxFile_t* file)
{
   polValidate(MTX_HERE, p);

   uint32_t hdr[3] = { MTX_TYPE_POLYNOMIAL, p->field, (uint32_t)p->degree };
   mfWrite32(file, hdr, 3);
   ffSetField(p->field);
   if (p->degree >= 0) {
      PTR tmp = ffAlloc(1, p->degree + 1);
      for (unsigned i = 0; i <= p->degree; ++i) {
         ffInsert(tmp, i, p->data[i]);
      }
      ffWriteRows(file, tmp, 1, p->degree + 1);
      ffFree(tmp);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes a single polynomial to a named file.
///
/// If a file with the same name alread eyists, its contents are destroyed.
/// See also @ref polWrite

void polSave(const Poly_t* pol, const char* fileName)
{
   polValidate(MTX_HERE, pol);
   MtxFile_t* file = mfOpen(fileName, "wb");
   polWrite(pol, file);
   mfClose(file);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads polynomial data from a file and return the polynomial.
///
/// Note: this function can only be called after a polynomial header was read. To simply read
/// a polynomial, use @ref polRead instead.

Poly_t *polReadData(MtxFile_t *f)
{
   const uint32_t objectType = mfObjectType(f);
   if (objectType != MTX_TYPE_POLYNOMIAL) {
      mtxAbort(MTX_HERE, "%s: bad type 0x%lx, expected 0x%lx (POLYNOMIAL)",
         f->name, (unsigned long) objectType, (unsigned long) MTX_TYPE_POLYNOMIAL);
   }

   const int field = (int) f->header[1];
   const size_t degree = f->header[2];
   ffSetField(field);
   Poly_t* polynomial = polAlloc(field, degree);
   if (polynomial->degree >= 0) {
      PTR tmpvec = ffAlloc(1, degree + 1);
      ffReadRows(f, tmpvec, 1, degree + 1);
      for (size_t i = 0; i <= degree; ++i) {
         polynomial->data[i] = ffExtract(tmpvec,i);
      }
      ffFree(tmpvec);
   }
   return polynomial;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads a polynomial from a file and returns the polynomial.
///
/// The function fails and aborts the program if no polynomial can be read. This include the case
/// that nothing could be read because the file pointer is already at the end of file.
/// To handle this case, use @ref mfTryReadHeader / @ref polReadData.

/// See also @ref polLoad.

Poly_t *polRead(MtxFile_t *f)
{
   mfReadHeader(f);
   return polReadData(f);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads a polynomial from a named file and returns the polynomial.
///
/// The first object in the file must be a polynomial, any further objects in the file are ignored.

Poly_t *polLoad(const char *fn)
{
   MtxFile_t* f = mfOpen(fn, "rb");
   Poly_t *p = polRead(f);
   mfClose(f);
   return p;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
