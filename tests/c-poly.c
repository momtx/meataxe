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
      p->data[i] = ffFromInt(mtxRandomInt(ffOrder));
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
      ASSERT(p[i]->degree == deg[i]);
      for (k = 0; k < p[i]->degree; ++k) {
         ASSERT(p[i]->data[k] == FF_ZERO);
      }
      if (p[i]->degree >= 0) {
         ASSERT(p[i]->data[p[i]->degree] == FF_ONE);
      }
   }
   for (i = 0; i < NPOLY; ++i) {
       polFree(p[i]);
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
      a->data[i] = FF_ONE;
      ASSERT(polCompare(a,b) != 0);
      b->data[i] = FF_ONE;
      ASSERT(polCompare(a,b) == 0);
   }
   if (ffGen != FF_ONE) {
      a->data[deg] = ffGen;
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
   ASSERT_EQ_INT(a->degree, 10);
   for (i = 0; i <= 10; ++i) {
      ASSERT_EQ_INT(a->data[i], FF_ONE);
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
   ASSERT_EQ_INT(a->degree,-1);
   polFree(a);
   polFree(b);

   /* Check x * 1 = x
      --------------- */
   a = polAlloc(ffOrder,1);
   b = polAlloc(ffOrder,0);
   polMul(a,b);
   ASSERT_EQ_INT(a->degree, 1);
   ASSERT_EQ_INT(a->data[0], FF_ZERO);
   ASSERT_EQ_INT(a->data[1], FF_ONE);
   polFree(a);
   polFree(b);

   /* Check (x+1)*(x^3-x) = x^4+x^3-x^2-x
      ----------------------------------- */
   a = polAlloc(ffOrder,1);
   b = polAlloc(ffOrder,3);
   a->data[0] = FF_ONE;
   b->data[1] = ffNeg(FF_ONE);
   polMul(a,b);
   ASSERT_EQ_INT(a->degree, 4);
   ASSERT_EQ_INT(a->data[0], FF_ZERO);
   ASSERT_EQ_INT(a->data[1], ffNeg(FF_ONE));
   ASSERT_EQ_INT(a->data[2], ffNeg(FF_ONE));
   ASSERT_EQ_INT(a->data[3], FF_ONE);
   ASSERT_EQ_INT(a->data[4], FF_ONE);
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
      ASSERT(a->degree <= 0 && b->degree <= 0);

      polFree(a);
      polFree(b);
      polFree(gcd);
      polFree(result[0]);
      polFree(result[1]);
      polFree(result[2]);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

struct FactorizationTestCase {
      int fieldOrder;
      int pol[6][20];   // degree, a[0], a[1], ... a[degree]
};


static int Polynomial_Factorization_(const struct FactorizationTestCase* tc, int mult)
{
   const size_t MAX_POLS = sizeof(tc->pol) / sizeof(tc->pol[0]);
   ffSetField(tc->fieldOrder);

   FPoly_t* expected = fpAlloc(tc->fieldOrder);
   Poly_t* poly = polAlloc(tc->fieldOrder, 0);
   for (size_t i = 0; i < MAX_POLS && tc->pol[i][0] > 0; ++i) {
      Poly_t* factor = polAlloc(tc->fieldOrder, tc->pol[i][0]);
      for (int k = 0; k <= factor->degree; ++k)
         factor->data[k] = ffFromInt(tc->pol[i][k+1]);
      for (int k = 0; k < mult; ++k) {
         polMul(poly, factor);
      }
      fpMulP(expected, factor, mult);
      polFree(factor);
   }

   FPoly_t* factorized = Factorization(poly);
   polFree(poly);
   ASSERT_EQ_INT(fpCompare(factorized, expected), 0);
   fpFree(factorized);
   fpFree(expected);
   return 0;
}


TstResult Polynomial_Factorization()
{
   const struct FactorizationTestCase TC[] = {
      { .fieldOrder = 243,
        .pol = {
           { 2, 2, 2, 1 },        // X2+2X+2
           { 3, 1, 2, 0, 1 },     // X3+2X+1
           { 4, 2, 0, 0, 2, 1 }   // X4+2X3+2
        } },
      { .fieldOrder = 256,
        .pol = {
           { 3, 1, 1, 0, 1 },             // X3+X+1
           { 5, 1, 0, 1, 0, 0, 1 },       // X5+X2+1
           { 7, 1, 1, 0, 0, 0, 0, 0, 1 }  // X7+X+1
        } },
   };
   const struct FactorizationTestCase* const TC_END = TC + sizeof(TC) / sizeof(TC[0]);

   int result = 0;
   for (const struct FactorizationTestCase* tc = TC; tc < TC_END; ++tc) {
      for (int mult = 1; mult <= 3; ++mult) {
         result |= Polynomial_Factorization_(tc, mult);
         if (result) { return result; }
      }
   }
   return result;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
