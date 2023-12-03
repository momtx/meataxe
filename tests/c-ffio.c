////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for ffio.c
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <string.h>

static const char* FILE_NAME = "test.tmp.1";

////////////////////////////////////////////////////////////////////////////////////////////////////

static void randomFill(PTR buf, size_t nor, size_t noc)
{
   for (size_t r = 0; r < nor; ++r) {
      PTR row = ffGetPtr(buf, r, noc);
      for (size_t c = 0; c < noc; ++c) {
         ffInsert(row, c, ffFromInt(mtxRandomInt(ffOrder)));
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int compareRows(PTR wrBuf, PTR rdBuf, size_t nor, size_t noc)
{
   ASSERT_EQ_INT(memcmp(rdBuf, wrBuf, nor * ffRowSize(noc)), 0);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void writeRowsWithBlockSize(
   MtxFile_t* file, PTR buf, uint32_t nor, uint32_t noc, uint32_t blockSize)
{
   for (uint32_t rowsWritten = 0; rowsWritten < nor;) {
      const uint32_t n = (rowsWritten + blockSize <= nor) ? blockSize : nor - rowsWritten;
      ffWriteRows(file, ffGetPtr(buf, rowsWritten, noc), n, noc);
      rowsWritten += n;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void readRowsWithBlockSize(
   MtxFile_t* file, PTR buf, uint32_t nor, uint32_t noc, uint32_t blockSize)
{
   for (uint32_t rowsRead = 0; rowsRead < nor;) {
      const uint32_t n = (rowsRead + blockSize <= nor) ? blockSize : nor - rowsRead;
      ffReadRows(file, ffGetPtr(buf, rowsRead, noc), n, noc);
      rowsRead += n;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_RowIo_CanWriteAndReadWithDifferentRowSizes()
{
   const uint32_t TBL_END = 0xFFFFFFFF;
   const uint32_t NOC[] = { 0, 1, 2, 31, 32, 33, 63, 64, 65, 127, 128, 129, TBL_END };
   const uint32_t BLOCK[] = { 1, 2, 3, 5, TBL_END };
   const size_t NROWS = 1000;

   int result = 0;
   MtxFile_t* file = mfOpen(FILE_NAME, "w+b");
   for (const uint32_t* pnoc = NOC; result == 0 && *pnoc != TBL_END; ++pnoc) {

      // Set up read and write buffers.
      PTR wrbuf = ffAlloc(NROWS, *pnoc);
      randomFill(wrbuf, NROWS, *pnoc);
      PTR rdbuf = ffAlloc(NROWS, *pnoc);

      for (const uint32_t* pblkw = BLOCK; result == 0 && *pblkw != TBL_END; ++pblkw) {
         sysFseek(file->file, 0);
         writeRowsWithBlockSize(file, wrbuf, NROWS, *pnoc, *pblkw);

         for (const uint32_t* pblkr = BLOCK; result == 0 && *pblkr != TBL_END; ++pblkr) {
            sysFseek(file->file, 0);
            readRowsWithBlockSize(file, rdbuf, NROWS, *pnoc, *pblkr);
            result |= compareRows(rdbuf, wrbuf, NROWS, *pnoc);
         }
      }

      ffFree(wrbuf);
      ffFree(rdbuf);
   }
   mfClose(file);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_RowIo_FailsOnPartialRead(int q)
{
   ffSetField(q);
   const size_t NOC = 10;
   const size_t NOR = 10;
   PTR rows = ffAlloc(NOR, NOC);

   // Write NOR - 1 rows
   MtxFile_t* file = mfOpen(FILE_NAME, "w+b");
   ffWriteRows(file, rows, NOR - 1, NOC);
   // Try to read NOR rows
   sysFseek(file->file, 0);
   ASSERT_ABORT(ffReadRows(file, rows, NOR, NOC));

   mfClose(file);
   remove(FILE_NAME);
   ffFree(rows);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestRowIo2(int noc)
{
   PTR rowBuf = ffAlloc(1, noc);

   MtxFile_t* file = mfOpen(FILE_NAME, "wb");
   ffWriteRows(file, rowBuf, 1, noc);
   for (int i = 0; i < noc; ++i) {
      ffInsert(rowBuf, i, FF_ONE);
      ffWriteRows(file, rowBuf, 1, noc);
   }
   mfClose(file);

   file = mfOpen(FILE_NAME, "rb");
   ffReadRows(file, rowBuf, 1, noc);
   for (int col = 0; col < noc; ++col) {
      ASSERT_EQ_INT(ffExtract(rowBuf, col), FF_ZERO);
   }
   for (int i = 0; i < noc; ++i) {
      ffReadRows(file, rowBuf, 1, noc);
      for (int col = 0; col <= i; ++col) {
         ASSERT_EQ_INT(ffExtract(rowBuf, col), FF_ONE);
      }
   }
   mfClose(file);

   remove(FILE_NAME);
   ffFree(rowBuf);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_RowIo(int q)
{
   int result = 0;
   for (int noc = 0; result == 0 && noc < 65; ++noc) {
      result |= TestRowIo2(noc);
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int CmpMat(PTR a, PTR b, int nor, int noc)
{
   while (nor-- > 0) {
      int diff = ffCmpRows(a,b, noc);
      if (diff != 0) {
         return diff;
      }
      ffStepPtr(&a, noc);
      ffStepPtr(&b, noc);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Test row i/o with header
static int TestHdr2(uint32_t nor, uint32_t noc, PTR buf1, PTR buf2)
{
   MtxFile_t* f;

   /* Write <buf1> into file
      ---------------------- */
   f = mfCreate(FILE_NAME,ffOrder,nor,noc);
   ffWriteRows(f,buf1,nor,noc);
   mfClose(f);

   /* Read the file header and check the value.
      ----------------------------------------- */
   memset(buf2,0,ffRowSize(noc) * nor);
   f = mfOpen(FILE_NAME, "rb");
   mfReadHeader(f);
   ASSERT_EQ_INT(f->header[0], ffOrder);
   ASSERT_EQ_INT(f->header[1], nor);
   ASSERT_EQ_INT(f->header[2], noc);

   // Read the rows. 
   ffReadRows(f,buf2,nor,noc);
   mfClose(f);

   // Compare <buf1> and <buf2>
   //tstPrintRows("buf1", buf1, nor, noc);
   //tstPrintRows("buf2", buf2, nor, noc);
   ASSERT_EQ_INT(CmpMat(buf1,buf2,nor,noc),0);

   remove(FILE_NAME);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_FileHeader(int q)
{
   int noc;
   const int nor = 100;

   int result = 0;
   for (noc = 0; noc < 65 && result == 0; ++noc) {
      PTR buf1, buf2;
      int i;
      uint32_t x = 0;

      buf1 = ffAlloc(nor, noc);
      buf2 = ffAlloc(nor, noc);
      for (i = 0; i < nor; ++i) {
         int k;
         PTR p = (PTR)((char *) buf1 + ffSize(i, noc));
         for (k = 0; k < noc; ++k) {
            ffInsert(p,k,ffFromInt((x >> 10) % ffOrder));
            x = x * 69069 + 13;
         }
      }

      result |= TestHdr2(nor,noc,buf1,buf2);

      free(buf1);
      free(buf2);
   }
   return result;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
