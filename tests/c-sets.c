#include "c-sets.h"
#include "check.h"
#include "meataxe.h"

#include <string.h>


MTX_DEFINE_FILE_INFO

static int ErrorFlag = 0;

static void MyErrorHandler(const MtxErrorRecord_t *err)
{
    ErrorFlag = 1;
    err = NULL;
}

static int CheckError()
{
    int i = ErrorFlag;
    ErrorFlag = 0;
    return i;
}


/* --------------------------------------------------------------------------
   TestSetAlloc() - Set allocation
   -------------------------------------------------------------------------- */

#define NMAT 5

void TestSetAlloc(unsigned flags)
{
    Set_t  *m[NMAT];
    MtxErrorHandler_t *old_err_handler;
    int i;

    for (i = 0; i < NMAT; ++i) 
	m[i] = SetAlloc();

    for (i = 0; i < NMAT; ++i) 
    {
	SetIsValid(m[i]);
	MTX_VERIFY(m[i]->Size == 0);
    }
    for (i = 0; i < NMAT; ++i) 
    {
	if (SetFree(m[i]) != 0) 
	    Error("SetFree() failed");
    }
    old_err_handler = MtxSetErrorHandler(MyErrorHandler);
    for (i = 0; i < NMAT; ++i) 
    {
	if (SetIsValid(m[i]) || !CheckError()) 
	    Error("SetIsValid() failed");
    }
    MtxSetErrorHandler(old_err_handler);
    flags = 0;
}



/* --------------------------------------------------------------------------
   TestSetOp() - Set insert/test
   -------------------------------------------------------------------------- */


void TestSetOp(unsigned flags)
{
    Set_t *s;
    long d[100];
    int i;


    MtxRandomInit(1213);
    memset(d,0,sizeof(d));
    for (i = 1; i <= 100; ++i)
    {
	int p;
	for (p = MtxRandomInt(100); d[p] != 0; p = (p + 1) % 100);
	d[p] = i;
    }

    s = SetAlloc();
    for (i = 0; i < 100; ++i)
    {
	int k;
	SetInsert(s,d[i]);
	if (s->Size != i + 1)
	    Error("Bad size");
	for (k = 0; k <= i; ++k)
	{
	    if (!SetContains(s,d[k]))
		Error("Element not inserted");
	}
	for (k = i + 1; k < 100; ++k)
	{
	    if (SetContains(s,d[k]))
		Error("Unexpected Element");
	}
    }
    SetFree(s);
    flags = 0;
}



