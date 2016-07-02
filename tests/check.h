////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests framework
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _CHECK_H_
#define _CHECK_H_

#include <stdarg.h>
#include <meataxe.h>

#if ZZZ == 0

#define ISFEL(f) ((unsigned int)(f) < (unsigned int)FfOrder)

#else

#define ISFEL(f) ((f) == 0xffff ||\
	(unsigned short)(f) < (unsigned short)FfOrder-1)
#error

#endif

typedef void test_F;
typedef test_F (*test_Function)();
typedef struct test_Definition {
   test_Function f;
   const char *name;
   const char *file;
   int line;
} test_Definition;


void test_Assert(const char *file, int line, const char *func, int e, const char *estr);
MTX_PRINTF_ATTRIBUTE(6,7)
void test_AssertF(const char *file, int line, const char *func, int e, const char *estr,
	const char *msg, ...);
#define ASSERT(e) test_Assert(__FILE__,__LINE__,__func__,e,#e)
#define ASSERT1(e,msg,a1) test_AssertF(__FILE__,__LINE__,__func__,e,#e,msg,a1)
#define ASSERT2(e,msg,a1,a2) test_AssertF(__FILE__,__LINE__,__func__,e,#e,msg,a1,a2)

void test_EqInt(const char *file, int line, const char *func, int act, int exp,
	        const char *actstr, const char *expstr);
#define ASSERT_EQ_INT(act,exp) test_EqInt(__FILE__,__LINE__,__func__,act,exp,#act,#exp)

MTX_PRINTF_ATTRIBUTE(4,5)
void test_Fail(const char *file, int line, const char *func, const char *msg, ...);
#define TST_FAIL(txt) test_Fail(__FILE__,__LINE__,__func__,"%s",txt)
#define TST_FAIL1(txt,a1) test_Fail(__FILE__,__LINE__,__func__,"%s",txt,a1)
#define TST_FAIL2(txt,a1,a2) test_Fail(__FILE__,__LINE__,__func__,"%s",txt,a1,a2)
#define TST_FAIL3(txt,a1,a2,a3) test_Fail(__FILE__,__LINE__,__func__,"%s",txt,a1,a2,a3)
#define TST_FAIL4(txt,a1,a2,a3,a4) test_Fail(__FILE__,__LINE__,__func__,"%s",txt,a1,a2,a3,a4)
#define TST_FAIL5(txt,a1,a2,a3,a4,a5) test_Fail(__FILE__,__LINE__,__func__,"%s",txt,a1,a2,a3,a4,a5)
#define TST_FAIL6(txt,a1,a2,a3,a4,a5,a6) \
   test_Fail(__FILE__,__LINE__,__func__,"%s",txt,a1,a2,a3,a4,a5,a6)

extern FEL *FTab;
extern int NextField();
extern void SelectField(int f);
extern void MakeFTab();
extern Matrix_t *MkMat(int nor, int noc, ...);

Perm_t *RndPerm(int degree);
Matrix_t *RndMat(int fl, int nor, int noc);
Poly_t *RndPol(int fl, int mindeg, int maxdeg);

void TstClearError();
void TstStartErrorChecking();
void TstStopErrorChecking();
int TstHasError();

#endif

