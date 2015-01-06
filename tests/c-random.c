#include "c-random.h"
#include "meataxe.h"
#include "check.h"

#include <string.h>





/* --------------------------------------------------------------------------
   Test1() - Check some random numbers
   -------------------------------------------------------------------------- */


static void Test11(unsigned seed, const long *table)
{
    int i;

    MtxRandomInit(seed);
    for (i = 0; i < 10; ++i, ++table)
    {
	int k;
	long val = MtxRandom() & 0x7FFFFFFF;
	if (val != *table)
	    Error("Got 0x%lx, expected 0x%lx",val,*table);
	for (k = 0; k < 61; ++k)
	    MtxRandom();
    }
}


static void Test1()
{
    static long Table0[10] = { 826837439, 1433481918, 1807203728, 1251143873, 
	498964889,886423565,167672701, 1728315981,248403305,1037767977 };
    static long Table1[10] = { 269167349, 1677366103,1597714250, 970290675, 
	436236141, 2108708678, 89648197, 1313827126, 514978688, 628812726 };
    Test11(0,Table0);
    Test11(1,Table1);
    Test11(0,Table0);
}



/* --------------------------------------------------------------------------
   Test2() - Check if the random numbers are sufficiently equally distributed
   -------------------------------------------------------------------------- */

static void Test21(int n, int *count)
{
    const int loop = 1500;
    const int min = loop - loop / 10;
    const int max = loop + loop / 10;

    int i;
    for (i = 0; i < n; ++i) count[i] = 0;
    for (i = 0; i < loop * n; ++i)
	++count[MtxRandomInt(n)];

    for (i = 0; i < n; ++i) 
    {
	if (count[i] < min || count[i] > max)
	    Error("Value %d hit %d times, expected %d..%d",i,count[i],min,max);
    }
}



static void Test2()
{
    int count[100];
    int n;
    for (n = 10; n < 100; ++n)
	Test21(n,count);
}




void TestRng(unsigned flags)
{
    Test1();
    Test2();
    flags = 0;
}

