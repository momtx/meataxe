////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Checks for various i/o functions.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult FileIo()
{
    Matrix_t *mat1, *mat2;
    Poly_t *pol1, *pol2;
    Perm_t *perm1, *perm2;
    FILE *f;

    SelectField(5);
    f = sysFopen("check.1","wb");
    mat1 = RndMat(5,30,30);
    matSave(mat1,"check.ma1");
    matWrite(mat1,f);
    pol1 = RndPol(5,100,200);
    polSave(pol1,"check.po1");
    polWrite(pol1,f);
    perm1 = RndPerm(100);
    permSave(perm1,"check.pe1");
    permWrite(perm1,f);
    fclose(f);

    mat2 = matLoad("check.ma1");
    ASSERT_EQ_INT(matCompare(mat1,mat2), 0);
    matFree(mat2);

    pol2 = polLoad("check.po1");
    ASSERT_EQ_INT(polCompare(pol1,pol2), 0);
    polFree(pol2);

    perm2 = permLoad("check.pe1");
    ASSERT_EQ_INT(permCompare(perm1,perm2), 0);
    permFree(perm2);

    f = sysFopen("check.1","rb");
    mat2 = matRead(f);
    pol2 = polRead(f);
    perm2 = permRead(f);
    fclose(f);

    ASSERT_EQ_INT(matCompare(mat1,mat2), 0);
    pol2 = polLoad("check.po1");
    ASSERT_EQ_INT(polCompare(pol1,pol2), 0);
    perm2 = permLoad("check.pe1");
    ASSERT_EQ_INT(permCompare(perm1,perm2), 0);

    matFree(mat1);
    polFree(pol1);
    permFree(perm1);
    matFree(mat2);
    polFree(pol2);
    permFree(perm2);
    return 0;
}
