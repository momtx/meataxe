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
//    Matrix_t *mat1, *mat2;
//    Poly_t *pol1, *pol2;
    Perm_t *perm1, *perm2;
    FILE *f;

    SelectField(5);
    f = sysFopen("check.1","wb");
//    mat1 = RndMat(5,30,30);
//    matSave(mat1,"check.ma1");
//    matWrite(mat1,f);
//    pol1 = RndPol(5,100,200);
//    polSave(pol1,"check.po1");
//    polWrite(pol1,f);
    perm1 = RndPerm(100);
    permSave(perm1,"check.pe1");
    permWrite(perm1,f);
    fclose(f);

//    mat2 = matLoad("check.ma1");
//    ASSERT_EQ_INT(matCompare(mat1,mat2), 0);
//    matFree(mat2);
//
//    pol2 = polLoad("check.po1");
//    ASSERT_EQ_INT(polCompare(pol1,pol2), 0);
//    polFree(pol2);

    perm2 = permLoad("check.pe1");
    ASSERT_EQ_INT(permCompare(perm1,perm2), 0);
    permFree(perm2);

//    f = sysFopen("check.1","rb");
//    mat2 = matRead(f);
//    pol2 = polRead(f);
//    perm2 = permRead(f);
//    fclose(f);
//
//    ASSERT_EQ_INT(matCompare(mat1,mat2), 0);
//    pol2 = polLoad("check.po1");
//    ASSERT_EQ_INT(polCompare(pol1,pol2), 0);
//    perm2 = permLoad("check.pe1");
//    ASSERT_EQ_INT(permCompare(perm1,perm2), 0);
//
//    matFree(mat1);
//    polFree(pol1);
//    permFree(perm1);
//    matFree(mat2);
//    polFree(pol2);
//    permFree(perm2);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static TstResult FileIo_Read16_(FILE *f)
{
    static const uint8_t DATA[] = {1, 2, 3, 4, 5, 6};
    fwrite(DATA, 1, sizeof(DATA), f);

    rewind(f);
    uint16_t data[3];
    sysRead16(f, data, 3);
    ASSERT_EQ_INT(data[0], 0x0201);
    ASSERT_EQ_INT(data[1], 0x0403);
    ASSERT_EQ_INT(data[2], 0x0605);
    return 0;

}

TstResult FileIo_Read16()
{
   static const char FILE_NAME[] = "test.data";
   FILE *f = sysFopen(FILE_NAME, "w+b");
   TstResult rv = FileIo_Read16_(f);
   fclose(f);
   remove(FILE_NAME);
   return rv;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static TstResult FileIo_Write16_(FILE *f)
{
    static const uint16_t DATA[3] = {0xaa55, 0xa55a, 0x3bb3};
    rewind(f);
    sysWrite16(f, DATA, 3);
    ASSERT(ftell(f) == 6);

    uint8_t data[6] = {0};
    rewind(f);
    ASSERT(fread(data, 1, sizeof(data), f) == sizeof(data));
    ASSERT(data[0] == 0x55);
    ASSERT(data[1] == 0xaa);
    ASSERT(data[2] == 0x5a);
    ASSERT(data[3] == 0xa5);
    ASSERT(data[4] == 0xb3);
    ASSERT(data[5] == 0x3b);
    return 0;

}

TstResult FileIo_Write16()
{
   static const char FILE_NAME[] = "test.data";
   FILE *f = sysFopen(FILE_NAME, "w+b");
   TstResult rv = FileIo_Write16_(f);
   fclose(f);
   remove(FILE_NAME);
   return rv;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static TstResult FileIo_Read32_(FILE *f)
{
    static const uint8_t DATA[] = {1, 2, 3, 4, 5, 6, 7, 8};
    fwrite(DATA, 1, sizeof(DATA), f);

    rewind(f);
    uint32_t data[2] = {0};
    sysRead32(f, data, 2);
    ASSERT_EQ_INT(data[0], 0x04030201);
    ASSERT_EQ_INT(data[1], 0x08070605);
    return 0;

}

TstResult FileIo_Read32()
{
   static const char FILE_NAME[] = "test.data";
   FILE *f = sysFopen(FILE_NAME, "w+b");
   TstResult rv = FileIo_Read32_(f);
   fclose(f);
   remove(FILE_NAME);
   return rv;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static TstResult FileIo_Write32_(FILE *f)
{
    static const uint32_t DATA[2] = {0x12345678, 0x87654321};
    rewind(f);
    sysWrite32(f, DATA, 2);
    ASSERT(ftell(f) == 8);

    uint8_t data[8];
    rewind(f);
    ASSERT(fread(data, 1, sizeof(data), f) == sizeof(data));
    ASSERT(data[0] == 0x78);
    ASSERT(data[1] == 0x56);
    ASSERT(data[2] == 0x34);
    ASSERT(data[3] == 0x12);
    ASSERT_EQ_INT(data[4],0x21);
    ASSERT_EQ_INT(data[5],0x43);
    ASSERT_EQ_INT(data[6],0x65);
    ASSERT_EQ_INT(data[7],0x87);
    return 0;
}

TstResult FileIo_Write32()
{
   static const char FILE_NAME[] = "test.data";
   FILE *f = sysFopen(FILE_NAME, "w+b");
   TstResult rv = FileIo_Write32_(f);
   fclose(f);
   remove(FILE_NAME);
   return rv;
}

