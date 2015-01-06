#ifndef CHECK_FUNCTION_TABLE

void TestSetAlloc(unsigned flags);
void TestSetOp(unsigned flags);

#else

  { 451, "Integer set allocation", TestSetAlloc },
  { 452, "Integer set operations", TestSetOp },
  
#endif
