////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check functions for matrices.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


////////////////////////////////////////////////////////////////////////////////////////////////////

Perm_t *RndPerm(int degree)
{
   int i;
   Perm_t *p = permAlloc(degree);
   for (i = 0; i < 2 * degree; ++i) {
      int a, b;
      a = mtxRandomInt(degree);
      b = mtxRandomInt(degree);
      if (a != b) {
         int c = p->data[a];
         p->data[a] = p->data[b];
         p->data[b] = c;
      }
   }
   return p;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Perm_AllocFree()
{
    #define NPERM 3
   static const int deg[NPERM] = { 0, 5, 700 };
   Perm_t *p[NPERM];
   int i;

   for (i = 0; i < NPERM; ++i) {
      p[i] = permAlloc(deg[i]);
   }
   for (i = 0; i < NPERM; ++i) {
      int k;
      permValidate(MTX_HERE, p[i]);
      ASSERT(p[i]->degree == deg[i]);
      for (k = 0; k < p[i]->degree; ++k) {
         ASSERT(p[i]->data[k] == k);
      }
   }
   for (i = 0; i < NPERM; ++i) {
      ASSERT(permFree(p[i]) == 0);
   }

   for (i = 0; i < NPERM; ++i) {
       ASSERT(!permIsValid(p[i]));
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static Perm_t *mkPerm(const int *a)
{
   int i;
   for (i = 0; a[i] >= 0; ++i);
   Perm_t *p = permAlloc(i);
   for (i = 0; i < p->degree; ++i)
      p->data[i] = a[i];
   return p;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int verifyOrder(const int *data, int order)
{
   Perm_t *p = mkPerm(data);
   //permPrint("perm", p);
   ASSERT_EQ_INT(permOrder(p), order);
   permFree(p);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Perm_Order()
{
   static const int p1[] = { 1,0,3,2,5,6,4, -1 };
   static const int p2[] = {-1};
   static const int p3[] = { 17,2,12,8,0,3,7,10,13,1,16,6,9,13,5,11,15,4,    -1 };
   static const int p4[] = { 0,2,3,4,1, -1 };

   int result = 0;
   result |= verifyOrder(p1,6);
   result |= verifyOrder(p2,1);
   result |= verifyOrder(p3,12);
   result |= verifyOrder(p4,4);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Perm_Multiply()
{
   static int p1[] = { 1,2,0,4,3, -1 };
   static int p2[] = { 0,1,3,2,4, -1 };
   static int p3[] = { 1,3,0,4,2, -1 };
   Perm_t *perm1, *perm2, *perm3;

   perm1 = mkPerm(p1);
   perm2 = mkPerm(p2);
   perm3 = mkPerm(p3);
   permMul(perm1,perm2);
   ASSERT(permCompare(perm1,perm3) == 0);
   permFree(perm1);
   permFree(perm2);
   permFree(perm3);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Perm_Power()
{
   Perm_t *p1;
   int i;

   p1 = RndPerm(1000);

   for (i = 0; i < 20; ++i) {
      Perm_t *p3 = permPower(p1,i);
      Perm_t *p2 = permAlloc(p1->degree);
      int k;
      for (k = 0; k < i; ++k) {
         permMul(p2,p1);
      }
      ASSERT(permCompare(p2,p3) == 0);
      permFree(p2);
      permFree(p3);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Perm_Inverse()
{
   int i;

   for (i = 0; i < 5000; ++i) {
      int k;
      Perm_t *p1 = RndPerm(i % 200 + 1);
      Perm_t *p2 = permInverse(p1);
      permMul(p2,p1);
      for (k = 0; k < p2->degree; ++k) {
         ASSERT(p2->data[k] == k);
      }
      permFree(p1);
      permFree(p2);
   }
   return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
