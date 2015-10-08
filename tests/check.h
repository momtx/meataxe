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
#define ASSERT(e) test_Assert(__FILE__,__LINE__,__func__,e,#e)
void test_EqInt(const char *file, int line, const char *func, int act, int exp,
	        const char *actstr, const char *expstr);
#define ASSERT_EQ_INT(act,exp) test_EqInt(__FILE__,__LINE__,__func__,act,exp,#act,#exp)


extern FEL *FTab;
extern void Error(char *msg, ...);
extern int NextField();
extern void SelectField(int f);
extern void MakeFTab();
extern Matrix_t *MkMat(int nor, int noc, ...);

Perm_t *RndPerm(int degree);
Matrix_t *RndMat(int fl, int nor, int noc);
Poly_t *RndPol(int fl, int mindeg, int maxdeg);


#endif

