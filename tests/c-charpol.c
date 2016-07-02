////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for characteristic polynomial
//
// (C) Copyright 1998-2016 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"

#include <stdlib.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static void CheckPoly(Poly_t *p, int degree, ...)
{
   int i;
   va_list al;
   va_start(al,degree);
   ASSERT_EQ_INT(p->Degree,degree);
   for (i = 0; i <= degree; ++i) {
      ASSERT_EQ_INT(p->Data[i],FfFromInt(va_arg(al,int)));
   }
   va_end(al);
   PolFree(p);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F CharacteristicPolynomial()
{
   Matrix_t *a;
   Poly_t *pf;

   SelectField(2);
   a = MkMat(6,6,
             1,0,0,0,0,0,    0,1,1,0,0,0,    0,0,0,0,1,0,
             0,0,1,1,0,0,    0,0,0,0,0,1,    0,0,0,0,1,1);

   pf = CharPolFactor(a);
   CheckPoly(pf,1,1,1);
   pf = CharPolFactor(NULL);
   CheckPoly(pf,4,0,1,0,0,1);
   pf = CharPolFactor(NULL);
   CheckPoly(pf,1,1,1);
   pf = CharPolFactor(NULL);
   ASSERT(pf == NULL);
   MatFree(a);
}
