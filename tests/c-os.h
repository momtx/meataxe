
#ifndef CHECK_FUNCTION_TABLE

void TestOs(unsigned flags);
void TestIntIo(unsigned flags);

#else

  { 221, "OS interface", TestOs },
  { 222, "File i/o: integer", TestIntIo },

#endif
