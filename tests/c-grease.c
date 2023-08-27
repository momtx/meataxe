////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check greased matrix operations.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestGrMapRow1(Matrix_t *m, int gr_level)
{
   Matrix_t *input = RndMat(ffOrder,m->Nor,m->Nor);
   GreasedMatrix_t *gm = GrMatAlloc(m,gr_level);
   PTR res_std = ffAlloc(1, m->Noc);
   PTR res_grease = ffAlloc(1, m->Noc);
   int i;

   for (i = 0; i < m->Nor; ++i) {
      PTR vec = matGetPtr(input,i);
      ffMapRow(vec,m->Data,m->Nor,m->Noc, res_std);
      GrMapRow(vec,gm,res_grease);
      ASSERT_EQ_INT(ffCmpRows(res_grease,res_std,m->Noc), 0);
   }
   sysFree(res_std);
   sysFree(res_grease);
   matFree(input);
   GrMatFree(gm);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult GreasedMapRow(int q)
{
      int result = 0;
#if MTX_ZZZ == 0
      int gr_level;
      int max_gr_level;
      int fpow;
      Matrix_t *m = RndMat(ffOrder,20,20);
      for (fpow = ffOrder, max_gr_level = 1;
            max_gr_level <= 16 && fpow < 66000;
            ++max_gr_level, fpow *= ffOrder) {
      }
      --max_gr_level;
      for (gr_level = 0; gr_level <= max_gr_level; ++gr_level) {
         result |= TestGrMapRow1(m,gr_level);
      }
      matFree(m);
#elif MTX_ZZZ == 1
   (void)TestGrMapRow1;
   printf("Greasing is not supported for ZZZ=1 - SKIPPING TEST\n");
#else
   #error
#endif
   return result;
}

// vim:sw=3:ts=8:et:cin
