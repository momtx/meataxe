#include "c-bitstring.h"
#include "check.h"
#include "meataxe.h"

#include <string.h>


MTX_DEFINE_FILE_INFO

static int ErrorFlag = 0;

static void MyErrorHandler(const MtxErrorRecord_t *err)
{
    err = NULL;
    ErrorFlag = 1;
}

#ifdef DEBUG

static int CheckError()
{
    int i = ErrorFlag;
    ErrorFlag = 0;
    return i;
}

#endif


/* --------------------------------------------------------------------------
   TestBsAlloc() - Matrix allocation
   -------------------------------------------------------------------------- */

#define NMAT 5

void TestBsAlloc(unsigned flags)
{
    static int bssize[NMAT] = { 0,1,10,100,1000 };
    BitString_t  *m[NMAT];
    MtxErrorHandler_t *old_err_handler;
    int i;

    flags = 0;	/* Just to avoid compiler warning */
    for (i = 0; i < NMAT; ++i) m[i] = BsAlloc(bssize[i]);
    for (i = 0; i < NMAT; ++i) 
    {
	int k;
	BsIsValid(m[i]);
	MTX_VERIFY(m[i]->Size == bssize[i]);
	for (k = 0; k < bssize[i]; ++k)
	{
	    if (BsTest(m[i],k) != 0)
		Error("New bit string not zero");
	}
    }
    for (i = 0; i < NMAT; ++i) 
	if (BsFree(m[i]) != 0) Error("BsFree() failed");
    old_err_handler = MtxSetErrorHandler(MyErrorHandler);
#ifdef DEBUG
    for (i = 0; i < NMAT; ++i) 
	if (BsIsValid(m[i]) || !CheckError()) Error("BsIsValid() failed");
#endif
    MtxSetErrorHandler(old_err_handler);
}



/* --------------------------------------------------------------------------
   TestBsOp() - Bit string operations
   -------------------------------------------------------------------------- */

static void TestSetClear(int size, BitString_t *a)

{
    int i;

    for (i = 0; i < size; ++i)
    {
	int k;
	for (k = 0; k < i; ++k)
	    if (!BsTest(a,k)) Error("Bit is still 0 after BsSet");
	for (; k < size; ++k)
	    if (BsTest(a,k)) Error("Unexpected set bit");
	BsSet(a,i);
    }
    for (i = 0; i < size; ++i)
    {
	int k;
	for (k = 0; k < i; ++k)
	    if (BsTest(a,k)) Error("Bit is still 1 after BsClear");
	for (; k < size; ++k)
	    if (!BsTest(a,k)) Error("Unexpected cleared bit");
	BsClear(a,i);
    }
   
    for (i = 0; i < size; ++i)
	BsSet(a,i);
    if (BsClearAll(a) != 0)
	Error("BsClearAll() failed");
    for (i = 0; i < size; ++i)
	if (BsTest(a,i)) Error("Bit %d is still 1 after BsClearAll",i);
}


void TestBsOp(unsigned flags)

{
    const int size = 50;
    BitString_t *a;

    flags = 0;	/* Avoid compiler warning */
    a = BsAlloc(size);
    TestSetClear(size,a);
    BsFree(a);
}


/* --------------------------------------------------------------------------
   TestCompare() - Bit string camparison
   -------------------------------------------------------------------------- */

static void TestCompare1(int size, BitString_t *a, BitString_t *b)

{
    int i;
    for (i = 0; i < size; ++i)
    {
	if (BsCompare(a,b) != 0) Error("BsCompare(x,x) != 0");
	BsSet(a,i);
	if (BsCompare(a,b) == 0) 
	    Error("BsCompare() did not detect difference in bit %d",i);
	BsSet(b,i);
    }
    for (i = 0; i < size; ++i)
    {
	if (BsCompare(a,b) != 0) Error("BsCompare(x,x) != 0");
	BsClear(a,i);
	if (BsCompare(a,b) == 0) 
	    Error("BsCompare() did not detect difference in bit %d",i);
	BsClear(b,i);
    }
}



void TestBsCompare(unsigned flags)

{
    const int size = 50;
    BitString_t *a, *b;
    a = BsAlloc(size);
    b = BsAlloc(size);
    TestCompare1(size,a,b);
    BsFree(b);
    BsFree(a);
    flags = 0;
}



/* --------------------------------------------------------------------------
   TestBsCopy() - Test BsCopy() and BsDup()
   -------------------------------------------------------------------------- */

static void TestCopy1(int size, BitString_t *a, BitString_t *b)

{
    int i;
    BitString_t *c;

    for (i = 0; i < size; i += 5)
	BsSet(a,i);
    BsCopy(b,a);
    if (BsCompare(a,b) != 0)
	Error("BsCopy() failed");
    c = BsDup(a);
    if (BsCompare(a,c) != 0)
	Error("BsDup() failed");
    BsFree(c);
}



void TestBsCopy(unsigned flags)
{
    const int size = 49;
    BitString_t *a, *b;
    a = BsAlloc(size);
    b = BsAlloc(size);
    TestCopy1(size,a,b);
    BsFree(b);
    BsFree(a);
    flags = 0;
}



/* --------------------------------------------------------------------------
   TestBsIo() - Test bit string file i/o
   -------------------------------------------------------------------------- */

static BitString_t *RndBs(int size)

{
    int i;
    BitString_t *bs = BsAlloc(size);
    for (i = 0; i < size; ++i)
    {
	if (MtxRandomInt(2) != 0)
	    BsSet(bs,i);
    }
    return bs;
}


static void CheckIo1(BitString_t **bs, int n)

{
    const char file_name[] = "check.1";
    FILE *f;
    int i;

    /* Write bit strings
       ----------------- */
    f = SysFopen(file_name,FM_CREATE);
    for (i = 0; i < n; ++i)
	BsWrite(bs[i],f);
    fclose(f);

    /* Read bit strings
       ---------------- */
    f = SysFopen(file_name,FM_READ);
    for (i = 0; i < n; ++i)
    {
	BitString_t *a = BsRead(f);
	if (BsCompare(a,bs[i]) != 0)
	    Error("Read error");
    }
    fclose(f);
}


void TestBsIo(unsigned flags)
{
    BitString_t *a[10];
    int i;

    MtxRandomInit(1235);
    a[0] = BsAlloc(0);
    for (i = 1; i < 10; ++i)
	a[i] = RndBs(MtxRandomInt(100));
    CheckIo1(a,10);
    for (i = 1; i < 10; ++i)
	BsFree(a[i]);
    flags = 0;
}



/* --------------------------------------------------------------------------
   TestBsAndOr() - Test bit string AND/OR/MINUS
   -------------------------------------------------------------------------- */

static void TestAndOr1(int size)

{
    BitString_t *a = BsAlloc(size);
    BitString_t *b = BsAlloc(size);
    BitString_t *bsor, *bsand, *bsminus;
    int i;

    for (i = 0; i < size; ++i)
    {
	if (i % 3 == 0) BsSet(a,i);
	if (i % 4 == 1) BsSet(b,i);
    }

    bsor = BsDup(a); BsOr(bsor,b);
    bsand = BsDup(a); BsAnd(bsand,b);
    bsminus = BsDup(a); BsMinus(bsminus,b);

    for (i = 0; i < size; ++i)
    {
	int exor = i % 3 == 0 || i % 4 == 1;
	int exand = i % 3 == 0 && i % 4 == 1;
	int exminus = i % 3 == 0 && i % 4 != 1;

	if (BsTest(bsor,i) ^ exor)
	    Error("Unexpected result in OR");
	if (BsTest(bsand,i) ^ exand)
	    Error("Unexpected result in AND");
	if (BsTest(bsminus,i) ^ exminus)
	    Error("Unexpected result in MINUS");
    }
    
    BsFree(bsor);
    BsFree(bsand);
    BsFree(bsminus);
    BsFree(b);
    BsFree(a);
}


void TestBsAndOr(unsigned flags)

{
    int i;

    for (i = 1; i < 100; i += i / 10 + 1)
	TestAndOr1(i);
    flags = 0;
}



/* --------------------------------------------------------------------------
   TestBsIntersectionCount() - Test BsIntersectionCount()
   -------------------------------------------------------------------------- */

static void CheckCount1(int size, BitString_t *a, BitString_t *b)

{
    int i;
    int count;
    int result;

    for (count = 0, i = 0; i < size; ++i)
    {
	if (BsTest(a,i) && BsTest(b,i))
	    ++count;
    }
    if ((result = BsIntersectionCount(a,b)) != count)
	Error("BsIntersectionCount()=%d, expected %d",result,count);
}


void TestBsIntersectionCount(unsigned flags)

{
    int i;
    MtxRandomInit(42);
    for (i = 0; i < 100; ++i)
    {
	int size = MtxRandomInt(100); 
	BitString_t *a = RndBs(size), *b = RndBs(size);
	CheckCount1(size,a,b);
	BsFree(a);
	BsFree(b);
    }
    flags = 0;
}




/* --------------------------------------------------------------------------
   TestBsIsSub() - Test BsIsSub()
   -------------------------------------------------------------------------- */

void TestBsIsSub(unsigned flags)

{
    const int size = 100;
    int i;
    BitString_t *a, *b;

    MtxRandomInit(42);
    for (i = 0; i < 100; ++i)
    {
	int k;
	a = RndBs(size);
	b = BsAlloc(size);
	if (!BsIsSub(b,a)) Error("Test 1 failed");
	BsCopy(b,a);
	if (!BsIsSub(b,a)) Error("Test 2 failed");
	for (k = 0; k < size; ++k)
	{
	    if (!BsTest(a,k))
	    {
		BsSet(b,k);
		if (BsIsSub(b,a)) Error("Test 3 failed");
		BsClear(b,k);
	    }
	}
	for (k = 0; k < size; ++k)
	{
	    BsClear(b,k);
	    if (!BsIsSub(b,a)) Error("Test 4 failed");
	}
        BsFree(a);
        BsFree(b);
    }
    flags = 0;
}


