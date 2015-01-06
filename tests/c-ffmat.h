#ifndef CHECK_FUNCTION_TABLE

void TestMapRow(unsigned flags);
void TestSumInter(unsigned flags);

#else

  { 141, "Vector-by matrix multiplication", TestMapRow },
  { 142, "Zasssenhaus algorithm", TestSumInter },

#endif
