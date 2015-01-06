#ifndef CHECK_FUNCTION_TABLE

void TestStf1(unsigned flags);
void TestStf2(unsigned flags);

#else

  { 750, "STF functions: basic", TestStf1 },
  { 751, "STF functions: Large data", TestStf2 },
  
#endif
