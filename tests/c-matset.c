/* ============================= C MeatAxe ==================================
   File:        $Id: c-matset.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Test functions for matrix sets.
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include "check.h"
#include "c-matrix.h"
#include "c-matset.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>



/* --------------------------------------------------------------------------
   TestMsClean() - Test MsClean
   -------------------------------------------------------------------------- */

static void TestMsClean1(MatrixSet_t *set)

{
    Matrix_t *zero;
    const int nor = 5;
    const int noc = 4;
    int i;

    for (i = 0; i < nor * noc; ++i)
    {
	int k;
	Matrix_t *m = MatAlloc(FfOrder,nor,noc);
	for (k = 0; k <= i; ++k)
	    FfInsert(MatGetPtr(m,k / noc),k % noc,FTab[k % (FfOrder-1) + 1]);
	if (MsCleanAndAppend(set,m) != 0)
	    Error("MsCleanAndAppend() failed");
    }

    zero = MatAlloc(FfOrder,nor,noc);
    for (i = 0; i < nor * noc; ++i)
    {
	Matrix_t *m = RndMat(FfOrder,nor,noc);
	MsClean(set,m);
	if (MatCompare(m,zero) != 0)
	    Error("MsClean() did not clean completely");
	MatFree(m);
    }
    MatFree(zero);
}




void TestMsClean(unsigned flags)

{
    MtxRandomInit(0);
    while (NextField() > 0)
    {
	MatrixSet_t *set;
	set = MsAlloc();
	TestMsClean1(set);
	MsFree(set);
    }
    flags = 0;
}








/* --------------------------------------------------------------------------
   TestMsAlloc() - Test matrix set creation/destruction
   -------------------------------------------------------------------------- */

void TestMsAlloc(unsigned flags)

{
    while (NextField() > 0)
    {
	Matrix_t *m1, *m2;
	MatrixSet_t *set;
	unsigned long magic;
	set = MsAlloc();
	if (set == NULL)
	    Error("MsAlloc() failed");
	m1 = RndMat(FfOrder,10,20);
	m2 = RndMat(FfOrder,10,20);
	magic = m1->Magic;
	if (   MsCleanAndAppend(set,m1) != 0
	    || MsCleanAndAppend(set,m2) != 0)
	    Error("MsCleanAndAppend() failed");
	if (MsFree(set) != 0)
	    Error ("MsFree() failed");
	if (m1->Magic == magic || m2->Magic == magic)
	    Error ("MsFree() did not destroy matrices");
    }
    flags = 0;
}




