#ifndef CHECK_FUNCTION_TABLE

void TestQuotProj(unsigned flags);
void TestQuotOp(unsigned flags);

#else

  { 801, "Quotient projection", TestQuotProj },
  { 802, "Quotient operation", TestQuotOp },

#endif
