#ifndef _CHECK_H_
#define _CHECK_H_

#include <stdarg.h>
#include "meataxe.h"

#if ZZZ == 0

#define ISFEL(f) ((unsigned int)(f) < (unsigned int)FfOrder)

#else

#define ISFEL(f) ((f) == 0xffff ||\
	(unsigned short)(f) < (unsigned short)FfOrder-1)
#error

#endif

extern FEL *FTab;
extern void Error(char *msg, ...);
extern int NextField();
extern void SelectField(int f);
extern void MakeFTab();
extern Matrix_t *MkMat(int nor, int noc, ...);

#endif

