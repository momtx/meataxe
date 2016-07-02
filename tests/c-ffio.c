////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for ffio.c
//
// (C) Copyright 1998-2016 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"

#include <stdlib.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestRowIo2(PTR row0, PTR row1, PTR buf)
{
   FILE *f;
   int i;
   long x;
   const char *file_name = "check.1";

   f = SysFopen(file_name,FM_CREATE);
   for (x = 0, i = 0; i < 100; ++i) {
      ASSERT_EQ_INT(FfWriteRows(f,(x & 0x1000) ? row0 : row1,1), 1);
      x = x * 69069 + 1;
   }
   fclose(f);

   f = SysFopen(file_name,FM_READ);
   for (x = 0, i = 0; i < 100; ++i) {
      ASSERT_EQ_INT(FfReadRows(f,buf,1), 1);
      if (FfCmpRows(buf,(x & 0x1000) ? row0 : row1) != 0) {
         TST_FAIL("Compare failed");
      }
      x = x * 69069 + 1;
   }
   fclose(f);

   remove(file_name);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestRowIo1()
{
   int noc;

   for (noc = 0; noc < 65; ++noc) {
      PTR row0, row1, buf;
      int i;

      FfSetNoc(noc);
      row0 = FfAlloc(1);
      row1 = FfAlloc(1);
      buf = FfAlloc(1);

      FfMulRow(row0,FF_ZERO);
      for (i = 0; i < noc - 1; ++i) {
         FfInsert(row1,i,FF_ONE);
      }
      TestRowIo2(row0,row1,buf);
      free(row0);
      free(row1);
      free(buf);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F RowIo()
{
   while (NextField() > 0) {
      TestRowIo1();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int CmpMat(PTR a, PTR b, int nor)
{
   while (nor-- > 0) {
      int diff = FfCmpRows(a,b);
      if (diff != 0) {
         return diff;
      }
      FfStepPtr(&a);
      FfStepPtr(&b);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Test row i/o with header
static void TestHdr2(int noc, PTR buf1, PTR buf2, int nor)
{
   FILE *f;
   const char *file_name = "check.1";
   int fld2 = -1, nor2 = -1, noc2 = -1;

   /* Write <buf1> into file
      ---------------------- */
   f = FfWriteHeader(file_name,FfOrder,nor,noc);
   ASSERT_EQ_INT(FfWriteRows(f,buf1,nor), nor);
   fclose(f);

   /* Read the file header and check the value.
      ----------------------------------------- */
   memset(buf2,0,FfRowSize(noc) * nor);
   f = FfReadHeader(file_name,&fld2,&nor2,&noc2);
   ASSERT_EQ_INT(fld2, FfOrder);
   ASSERT_EQ_INT(nor2, nor);
   ASSERT_EQ_INT(noc2, noc);

   // Read the rows. If <noc> is not zero, we try to read one more row to check if
   // <FfreadRows()> handles the EOF correctly. For <noc> = 0 FfReadRows() always returns
   // the requested number of rows, so the check is not possible.
   ASSERT_EQ_INT(FfReadRows(f,buf2,(noc == 0) ? nor : nor + 1), nor);
   fclose(f);

   // Compare <buf1> and <buf2>
   ASSERT_EQ_INT(CmpMat(buf1,buf2,nor),0);

   remove(file_name);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestHdr1()
{
   int noc;
   const int bufsize = 100;

   for (noc = 0; noc < 65; ++noc) {
      PTR buf1, buf2;
      int i;
      long x = 0;

      FfSetNoc(noc);
      buf1 = FfAlloc(bufsize);
      buf2 = FfAlloc(bufsize);
      for (i = 0; i < bufsize; ++i) {
         int k;
         PTR p = (PTR)((char *) buf1 + i * FfCurrentRowSize);
         for (k = 0; k < noc; ++k) {
            FfInsert(p,k,FfFromInt(x >> 10));
            x = x * 69069 + 13;
         }
      }

      TestHdr2(noc,buf1,buf2,bufsize);

      free(buf1);
      free(buf2);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F FileHeader()
{
   while (NextField() > 0) {
      TestHdr1();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestSeek2(int noc, PTR buf1, PTR buf2, int nor)
{
   FILE *f;
   const char *file_name = "check.1";
   int fld2, nor2, noc2;
   int i;

   /* Write <buf1> into file
      ---------------------- */
   f = FfWriteHeader(file_name,FfOrder,nor,noc);
   ASSERT_EQ_INT(FfWriteRows(f,buf1,nor), nor);
   fclose(f);

   /* Read the rows in reverse order.
      ------------------------------- */
   memset(buf2,0,FfRowSize(noc) * nor);
   f = FfReadHeader(file_name,&fld2,&nor2,&noc2);
   for (i = nor - 1; i >= 0; --i) {
      PTR x = (PTR)((char *)buf2 + i * FfCurrentRowSize);
      FfSeekRow(f,i);
      ASSERT_EQ_INT(FfReadRows(f,x,1), 1);
   }
   fclose(f);

   /* Compare <buf1> and <buf2>
      ------------------------- */
   ASSERT_EQ_INT(CmpMat(buf1,buf2,nor), 0);

   remove(file_name);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestSeek1()
{
   int noc;
   const int bufsize = 100;

   for (noc = 0; noc < 65; ++noc) {
      PTR buf1, buf2;
      int i;
      long x = 0;

      FfSetNoc(noc);
      buf1 = FfAlloc(bufsize);
      buf2 = FfAlloc(bufsize);
      for (i = 0; i < bufsize; ++i) {
         int k;
         PTR p = (PTR)((char *)buf1 + i * FfCurrentRowSize);
         for (k = 0; k < noc; ++k) {
            FfInsert(p,k,FfFromInt(x >> 10));
            x = x * 69069 + 13;
         }
      }

      TestSeek2(noc,buf1,buf2,bufsize);

      free(buf1);
      free(buf2);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F Seek()
{
   while (NextField() > 0) {
      TestSeek1();
   }
}
