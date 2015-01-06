#ifndef CHECK_FUNCTION_TABLE

void TestFelToInt(unsigned flags);
void TestField(unsigned flags);
void TestGen(unsigned flags);
void TestInsertExtract(unsigned flags);
void TestFindPiv(unsigned flags);
void TestSubfields(unsigned flags);
void TestRowOps(unsigned flags);
void TestPtr(unsigned flags);
void TestRowSize(unsigned flags);
void TestMulRow(unsigned flags);
void TestCmpRows(unsigned flags);

#else

  { 111, "Row operations", TestRowOps },
  { 101, "Compare operations", TestCmpRows },
  { 102, "FEL <--> integer conversion", TestFelToInt },
  { 103, "Finite field arithmetic", TestField },
  { 104, "Finite field generator", TestGen },
  { 105, "Subfield embedding/restriction", TestSubfields },
  { 106, "Row sizes", TestRowSize },
  { 107, "Row pointer arithmetic", TestPtr },
  { 108, "Insert/Extract", TestInsertExtract },
  { 109, "Finding pivot elements", TestFindPiv },
  { 110, "Row by scalar multiplication", TestMulRow },

#endif
