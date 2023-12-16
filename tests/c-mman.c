////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Textx for the memory management (mm) module.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Mm_CanRollback_NoObjects()
{
    uint32_t checkpoint = mmCheckpoint();
    permFree(permAlloc(100));
    mmRollback(checkpoint);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Mm_CanRollback_MultipleObjects()
{
    Perm_t* p0 = permAlloc(10);
    uint32_t checkpoint = mmCheckpoint();
    Perm_t* p1 = permAlloc(20);
    Poly_t* pol1 = polAlloc(2,10);
    Matrix_t* m1 = matAlloc(2,11,11);
    BitString_t* bs1 = bsAllocEmpty();
    mmRollback(checkpoint);
    ASSERT(!permIsValid(p1));
    ASSERT(!bsIsValid(bs1));
    ASSERT(!matIsValid(m1));
    ASSERT(!polIsValid(pol1));
    permFree(p0);
    return 0;
}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
