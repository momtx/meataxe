#ifndef CHECK_FUNCTION_TABLE

void TestRowIo(unsigned flags);
void TestHdr(unsigned flags);
void TestSeek(unsigned flags);

#else

  { 131, "FfSeekRow()", TestSeek },
  { 132, "Finite field file i/o", TestRowIo },
  { 133, "MeatAxe file format", TestHdr },

#endif
