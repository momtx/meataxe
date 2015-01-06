#ifndef CHECK_FUNCTION_TABLE

void TestBsAlloc(unsigned flags);
void TestBsAndOr(unsigned flags);
void TestBsCompare(unsigned flags);
void TestBsCopy(unsigned flags);
void TestBsOp(unsigned flags);
void TestBsIntersectionCount(unsigned flags);
void TestBsIo(unsigned flags);
void TestBsIsSub(unsigned flags);

#else

  { 401, "Bit string allocation", TestBsAlloc },
  { 402, "Bit manipulations", TestBsOp },
  { 403, "Bit string comparison", TestBsCompare },
  { 404, "Bit string copying", TestBsCopy },
  { 405, "Bit string file i/o", TestBsIo },
  { 406, "Bit string intersection count", TestBsIntersectionCount },
  { 407, "Bit string operations", TestBsAndOr },
  { 408, "Bit string incidence", TestBsIsSub },

#endif
