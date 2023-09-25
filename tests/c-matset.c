////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Test functions for matrix sets.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult MatrixSet_Clean(int q)
{
    MatrixSet_t *set;
    set = msAlloc();

    Matrix_t *zero;
    const int nor = 5;
    const int noc = 4;
    int i;

    for (i = 0; i < nor * noc; ++i) {
	int k;
	Matrix_t *m = matAlloc(ffOrder,nor,noc);
	for (k = 0; k <= i; ++k) {
	    ffInsert(matGetPtr(m,k / noc),k % noc,FTab[k % (ffOrder - 1) + 1]);
	}
	ASSERT_EQ_INT(msCleanAndAppend(set,m), 0);
    }

    zero = matAlloc(ffOrder,nor,noc);
    for (i = 0; i < nor * noc; ++i) {
	Matrix_t *m = RndMat(ffOrder,nor,noc);
	msClean(set,m);
	ASSERT_EQ_INT(matCompare(m,zero), 0);
	matFree(m);
    }
    matFree(zero);

    msFree(set);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult MatrixSet_Allocation(int q)
{
    Matrix_t *m1, *m2;
    unsigned long magic;
    MatrixSet_t *set = msAlloc();
    ASSERT(set != NULL);
    m1 = RndMat(ffOrder,10,20);
    m2 = RndMat(ffOrder,10,20);
    magic = m1->typeId;
    ASSERT(msCleanAndAppend(set,m1) == 0);
    ASSERT(msCleanAndAppend(set,m2) == 0);
    ASSERT(msFree(set) == 0);
    ASSERT(m1->typeId != magic);
    ASSERT(m2->typeId != magic);
    return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
