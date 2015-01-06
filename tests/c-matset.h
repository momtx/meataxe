
#ifndef CHECK_FUNCTION_TABLE

#include "meataxe.h"

void TestMsAlloc(unsigned flags);
void TestMsClean(unsigned flags);


#else

  { 340, "Matrix set creation/destruction", TestMsAlloc },
  { 341, "Matrix set cleaning", TestMsClean },

#endif
