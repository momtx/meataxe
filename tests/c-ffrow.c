////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check low-level row operations.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"

#include <stdlib.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestScalarProduct1(PTR a, PTR b, int size)
{
   int count;
   for (count = 0; count < 10; ++count) {
      int i;
      FEL expected = FF_ZERO, found;
      FfMulRow(a,FF_ZERO);
      FfMulRow(b,FF_ZERO);
      for (i = 0; i < size; ++i) {
         FEL f1 = FTab[MtxRandomInt(FfOrder)];
         FEL f2 = FTab[MtxRandomInt(FfOrder)];
         FfInsert(a,i,f1);
         FfInsert(b,i,f2);
         expected = FfAdd(expected,FfMul(f1,f2));
      }
      found = FfScalarProduct(a,b);
      if (found != expected) {
         Error("FfScalarProduct() returned %d, expected %d",found,expected);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
 
test_F ScalarProduct()
{
   MtxRandomInit(1231);
   while (NextField() > 0) {
      int size;
      for (size = 0; size < 1000; size += size / 10 + 1) {
         PTR a, b;
         FfSetNoc(size);
         a = FfAlloc(1);
         b = FfAlloc(1);
         TestScalarProduct1(a,b,size);
         SysFree(a);
         SysFree(b);
      }
   }
}
