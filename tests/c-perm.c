////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check functions for matrices.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

MTX_DEFINE_FILE_INFO

static int ErrorFlag = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////

static void MyErrorHandler(const MtxErrorRecord_t *err)
{
   ErrorFlag = 1;
   err = NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int CheckError()
{
   int i = ErrorFlag;
   ErrorFlag = 0;
   return i;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Perm_t *RndPerm(int degree)
{
   int i;
   Perm_t *p = PermAlloc(degree);
   for (i = 0; i < 2 * degree; ++i) {
      int a, b;
      a = MtxRandomInt(degree);
      b = MtxRandomInt(degree);
      if (a != b) {
         int c = p->Data[a];
         p->Data[a] = p->Data[b];
         p->Data[b] = c;
      }
   }
   return p;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

#define NPERM 3

test_F PermutationAlloc()
{
   static const int deg[NPERM] = { 0, 5, 700 };
   Perm_t *p[NPERM];
   MtxErrorHandler_t *old_err_handler;
   int i;

   for (i = 0; i < NPERM; ++i) {
      p[i] = PermAlloc(deg[i]);
   }
   for (i = 0; i < NPERM; ++i) {
      int k;
      PermIsValid(p[i]);
      MTX_VERIFY(p[i]->Degree == deg[i]);
      for (k = 0; k < p[i]->Degree; ++k) {
         MTX_VERIFY(p[i]->Data[k] == k);
      }
   }
   for (i = 0; i < NPERM; ++i) {
      if (PermFree(p[i]) != 0) {
	  TST_FAIL("PermFree() failed");
      }
   }
   old_err_handler = MtxSetErrorHandler(MyErrorHandler);
   for (i = 0; i < NPERM; ++i) {
      if (PermIsValid(p[i]) || !CheckError()) {
	  TST_FAIL("PermIsValid() failed");
      }
   }
   MtxSetErrorHandler(old_err_handler);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestOrd(long *data, int order)
{
   int deg;
   Perm_t *p;

   for (deg = 0; data[deg] >= 0; ++deg);
   p = PermAlloc(deg);
   memcpy(p->Data,data,sizeof(long) * deg);
   ASSERT_EQ_INT(PermOrder(p), order);
   PermFree(p);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F PermutationOrder()
{
   long p1[] = { 1,0,3,2,5,6,4, -1 };
   long p2[] = {-1};
   long p3[] = { 17,2,12,8,0,3,7,10,13,1,16,6,9,13,5,11,15,4,    -1 };
   long p4[] = { 0,2,3,4,1, -1 };

   TestOrd(p1,6);
   TestOrd(p2,1);
   TestOrd(p3,12);
   TestOrd(p4,4);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static Perm_t *MkPerm(int *a)
{
   int i;
   Perm_t *p;

   for (i = 0; a[i] >= 0; ++i) {
   }
   p = PermAlloc(i);
   for (i = 0; i < p->Degree; ++i) {
      p->Data[i] = a[i];
   }
   return p;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F PermutationMultiply()
{
   static int p1[] = { 1,2,0,4,3, -1 };
   static int p2[] = { 0,1,3,2,4, -1 };
   static int p3[] = { 1,3,0,4,2, -1 };
   Perm_t *perm1, *perm2, *perm3;

   perm1 = MkPerm(p1);
   perm2 = MkPerm(p2);
   perm3 = MkPerm(p3);
   PermMul(perm1,perm2);
   if (PermCompare(perm1,perm3) != 0) {
      TST_FAIL("PermMul() failed");
   }
   PermFree(perm1);
   PermFree(perm2);
   PermFree(perm3);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F PermutationPower()
{
   Perm_t *p1;
   int i;

   p1 = RndPerm(1000);

   for (i = 0; i < 20; ++i) {
      Perm_t *p3 = PermPower(p1,i);
      Perm_t *p2 = PermAlloc(p1->Degree);
      int k;
      for (k = 0; k < i; ++k) {
         PermMul(p2,p1);
      }
      if (PermCompare(p2,p3) != 0) {
         TST_FAIL("PermPwr() failed");
      }
      PermFree(p2);
      PermFree(p3);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F PermutationInverse()
{
   int i;

   for (i = 0; i < 5000; ++i) {
      int k;
      Perm_t *p1 = RndPerm(i % 200 + 1);
      Perm_t *p2 = PermInverse(p1);
      PermMul(p2,p1);
      for (k = 0; k < p2->Degree; ++k) {
         if (p2->Data[k] != k) {
            TST_FAIL("Inverse of permutation failed");
         }
      }
      PermFree(p1);
      PermFree(p2);
   }
}
