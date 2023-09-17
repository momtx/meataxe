////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check functions for matrices.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


////////////////////////////////////////////////////////////////////////////////////////////////////

Poly_t *RndPol(int fl, int mindeg, int maxdeg)
{
   Poly_t *p;
   int deg;
   int i;

   deg = mtxRandomInt(maxdeg - mindeg + 1) + mindeg;
   p = polAlloc(fl,deg);
   for (i = 0; i <= deg - 1; ++i) {
      p->Data[i] = ffFromInt(mtxRandomInt(ffOrder));
   }
   return p;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

#define NPOLY 5

TstResult Polynomial_Alloc(int q)
{
   static const int deg[NPOLY] = { -1,0,5,10,200 };
   Poly_t *p[NPOLY];
   int i;

   for (i = 0; i < NPOLY; ++i) {
      p[i] = polAlloc(q,deg[i]);
   }
   for (i = 0; i < NPOLY; ++i) {
      int k;
      ASSERT(polIsValid(p[i]));
      ASSERT(p[i]->Degree == deg[i]);
      for (k = 0; k < p[i]->Degree; ++k) {
         ASSERT(p[i]->Data[k] == FF_ZERO);
      }
      if (p[i]->Degree >= 0) {
         ASSERT(p[i]->Data[p[i]->Degree] == FF_ONE);
      }
   }
   for (i = 0; i < NPOLY; ++i) {
       ASSERT(polFree(p[i]) == 0);
   }
      
   for (i = 0; i < NPOLY; ++i) {
       ASSERT(!polIsValid(p[i]));
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Polynomial_AbortsOnDoubleFree()
{
   Poly_t *pol = polAlloc(3, 10);
   ASSERT(polIsValid(pol));
   polFree(pol);
   ASSERT(!polIsValid(pol));
   ASSERT_ABORT(polFree(pol));
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

#define NPOLY 5

static int TestPolCompare2(int fl1, int deg1, int fl2, int deg2, int result)
{
   Poly_t *a, *b;
   a = polAlloc(fl1,deg1);
   b = polAlloc(fl2,deg2);
   ASSERT_EQ_INT(polCompare(a,b),result);
   ASSERT_EQ_INT(polCompare(b,a),-result);
   polFree(a);
   polFree(b);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Polynomial_Compare2()
{
   int result = 0;
   result |= TestPolCompare2(2,-1,3,-1,-1);
   result |= TestPolCompare2(2, 0,3,-1,-1);
   result |= TestPolCompare2(2,10,3, 0,-1);
   result |= TestPolCompare2(2,0,2,-1,1);
   result |= TestPolCompare2(2,10,2,9,1);
   result |= TestPolCompare2(3,0,3,0,0);
   result |= TestPolCompare2(3,-1,3,-1,0);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestPolCompare1(int fl)
{
   Poly_t *a, *b;
   const int deg = 10;
   int i;

   a = polAlloc(fl,deg);
   b = polAlloc(fl,deg);
   ASSERT(polCompare(a,b) == 0);
   for (i = 0; i < deg; ++i) {
      a->Data[i] = FF_ONE;
      ASSERT(polCompare(a,b) != 0);
      b->Data[i] = FF_ONE;
      ASSERT(polCompare(a,b) == 0);
   }
   if (ffGen != FF_ONE) {
      a->Data[deg] = ffGen;
      ASSERT(polCompare(a,b) != 0);
   }
   polFree(a);
   polFree(b);
   return 0;
}

TstResult Polynomial_Compare1(int q)
{
   return TestPolCompare1(ffOrder);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestPolAdd1(int fl)
{
   Poly_t *a, *b, *c;

   a = polAlloc(fl,-1);
   b = polAlloc(fl,0);
   c = polAlloc(fl,0);
   polAdd(a,b);
   ASSERT(polCompare(a,c) == 0);
   polFree(a);
   polFree(b);
   polFree(c);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestPolAdd2(int fl)
{
   Poly_t *a, *b;
   int i;

   a = polAlloc(fl,-1);
   for (i = -1; i <= 10; ++i) {
      b = polAlloc(fl,i);
      polAdd(a,b);
      polFree(b);
   }
   ASSERT_EQ_INT(a->Degree, 10);
   for (i = 0; i <= 10; ++i) {
      ASSERT_EQ_INT(a->Data[i], FF_ONE);
   }
   polFree(a);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult PolynomialAdd(int q)
{
   int result = 0;
   result |= TestPolAdd1(ffOrder);
   result |= TestPolAdd2(ffOrder);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestPolMul1()
{
   Poly_t *a, *b;

   /* Check x * 0 = 0
      --------------- */
   a = polAlloc(ffOrder,1);
   b = polAlloc(ffOrder,-1);
   polMul(a,b);
   ASSERT_EQ_INT(a->Degree,-1);
   polFree(a);
   polFree(b);

   /* Check x * 1 = x
      --------------- */
   a = polAlloc(ffOrder,1);
   b = polAlloc(ffOrder,0);
   polMul(a,b);
   ASSERT_EQ_INT(a->Degree, 1);
   ASSERT_EQ_INT(a->Data[0], FF_ZERO);
   ASSERT_EQ_INT(a->Data[1], FF_ONE);
   polFree(a);
   polFree(b);

   /* Check (x+1)*(x^3-x) = x^4+x^3-x^2-x
      ----------------------------------- */
   a = polAlloc(ffOrder,1);
   b = polAlloc(ffOrder,3);
   a->Data[0] = FF_ONE;
   b->Data[1] = ffNeg(FF_ONE);
   polMul(a,b);
   ASSERT_EQ_INT(a->Degree, 4);
   ASSERT_EQ_INT(a->Data[0], FF_ZERO);
   ASSERT_EQ_INT(a->Data[1], ffNeg(FF_ONE));
   ASSERT_EQ_INT(a->Data[2], ffNeg(FF_ONE));
   ASSERT_EQ_INT(a->Data[3], FF_ONE);
   ASSERT_EQ_INT(a->Data[4], FF_ONE);
   polFree(a);
   polFree(b);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestPolMul2()
{
   int i;
   for (i = 0; i < 20; ++i) {
      Poly_t *a, *b, *c, *ab, *ba;
      a = RndPol(ffOrder,0,100);
      b = RndPol(ffOrder,0,100);
      c = RndPol(ffOrder,0,100);

      ab = polDup(a);
      polMul(ab,b);
      ba = polDup(b);
      polMul(ba,a);
      ASSERT_EQ_INT(polCompare(ab,ba), 0);

      polMul(ab,c);
      polMul(b,c);
      polMul(a,b);
      ASSERT_EQ_INT(polCompare(a,ab), 0);

      polFree(a);
      polFree(b);
      polFree(c);
      polFree(ab);
      polFree(ba);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult PolynomialMultiply(int q)
{
   int result = 0;
   result |= TestPolMul1();
   result |= TestPolMul2();
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Polynomial_Gcd()
{
   int n;

   for (n = 100; n > 0; --n) {
      Poly_t *a, *b, *gcd, *result[3];

      /* Create two random polynomials
         ----------------------------- */
      a = RndPol(ffOrder,0,100);
      b = RndPol(ffOrder,0,100);

      /* Calculate the g.c.d.
         -------------------- */
      ASSERT((gcd = polGcd(a,b)) != NULL);
      ASSERT(polGcdEx(a,b,result) == 0);

      /* Check the result
         ---------------- */
      ASSERT(polCompare(gcd,result[0]) == 0);
      
      polMul(result[1],a);
      polMul(result[2],b);
      polAdd(result[1],result[2]);
      ASSERT(polCompare(gcd,result[1]) == 0);
      ASSERT(polMod(a,gcd) != NULL);
      ASSERT(polMod(b,gcd) != NULL);
      ASSERT(a->Degree <= 0 && b->Degree <= 0);

      polFree(a);
      polFree(b);
      polFree(gcd);
      polFree(result[0]);
      polFree(result[1]);
      polFree(result[2]);
   }
   return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
