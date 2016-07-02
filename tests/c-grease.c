////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check greased matrix operations.
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

static void TestGrMapRow1(Matrix_t *m, int gr_level)
{
   Matrix_t *input = RndMat(FfOrder,m->Nor,m->Nor);
   GreasedMatrix_t *gm = GrMatAlloc(m,gr_level);
   PTR res_std = FfAlloc(1);
   PTR res_grease = FfAlloc(1);
   int i;

   for (i = 0; i < m->Nor; ++i) {
      PTR vec = MatGetPtr(input,i);
      FfSetNoc(m->Noc);
      FfMapRow(vec,m->Data,m->Nor,res_std);
      GrMapRow(vec,gm,res_grease);
      ASSERT_EQ_INT(FfCmpRows(res_grease,res_std), 0);
   }
   SysFree(res_std);
   SysFree(res_grease);
   MatFree(input);
   GrMatFree(gm);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F GreasedMapRow()
{
   while (NextField() > 0) {
      int gr_level;
      int max_gr_level;
      int fpow;
      Matrix_t *m = RndMat(FfOrder,20,20);
      for (fpow = FfOrder, max_gr_level = 1;
           max_gr_level <= 16 && fpow < 66000;
           ++max_gr_level, fpow *= FfOrder) {
      }
      --max_gr_level;
      for (gr_level = 0; gr_level <= max_gr_level; ++gr_level) {
         TestGrMapRow1(m,gr_level);
      }
      MatFree(m);
   }
}
