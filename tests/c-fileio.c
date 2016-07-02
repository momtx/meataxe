////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Checks for various i/o functions.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"

#include <stdlib.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F FileIo()
{
   Matrix_t *mat1, *mat2;
   Poly_t *pol1, *pol2;
   Perm_t *perm1, *perm2;
   FILE *f;

   SelectField(5);
   f = SysFopen("check.1",FM_CREATE);
   mat1 = RndMat(5,30,30);
   MatSave(mat1,"check.ma1");
   MatWrite(mat1,f);
   pol1 = RndPol(5,100,200);
   PolSave(pol1,"check.po1");
   PolWrite(pol1,f);
   perm1 = RndPerm(100);
   PermSave(perm1,"check.pe1");
   PermWrite(perm1,f);
   fclose(f);

   mat2 = MatLoad("check.ma1");
   ASSERT_EQ_INT(MatCompare(mat1,mat2), 0);
   MatFree(mat2);

   pol2 = PolLoad("check.po1");
   ASSERT_EQ_INT(PolCompare(pol1,pol2), 0);
   PolFree(pol2);

   perm2 = PermLoad("check.pe1");
   ASSERT_EQ_INT(PermCompare(perm1,perm2), 0);
   PermFree(perm2);

   f = SysFopen("check.1",FM_READ);
   mat2 = MatRead(f);
   pol2 = PolRead(f);
   perm2 = PermRead(f);
   fclose(f);

   ASSERT_EQ_INT(MatCompare(mat1,mat2), 0);
   pol2 = PolLoad("check.po1");
   ASSERT_EQ_INT(PolCompare(pol1,pol2), 0);
   perm2 = PermLoad("check.pe1");
   ASSERT_EQ_INT(PermCompare(perm1,perm2), 0);

   MatFree(mat1);
   PolFree(pol1);
   PermFree(perm1);
   MatFree(mat2);
   PolFree(pol2);
   PermFree(perm2);
}
