////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for ffio.c
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestRowIo2(PTR row0, PTR row1, PTR buf, int noc)
{
   FILE *f;
   int i;
   long x;
   const char *file_name = "check.1";

   f = sysFopen(file_name,"wb");
   for (x = 0, i = 0; i < 100; ++i) {
      ffWriteRows(f,(x & 0x1000) ? row0 : row1, 1, noc);
      x = x * 69069 + 1;
   }
   fclose(f);

   f = sysFopen(file_name,"rb");
   for (x = 0, i = 0; i < 100; ++i) {
      ffReadRows(f,buf,1,noc);
      ASSERT(ffCmpRows(buf,(x & 0x1000) ? row0 : row1, noc) == 0);
      x = x * 69069 + 1;
   }
   fclose(f);

   remove(file_name);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_RowIo(int q)
{
   int result = 0;
   for (int noc = 0; result == 0 && noc < 65; ++noc) {

      PTR row0, row1, buf;
      int i;

      row0 = ffAlloc(1, noc);
      row1 = ffAlloc(1, noc);
      buf = ffAlloc(1, noc);

      ffMulRow(row0,FF_ZERO, noc);
      for (i = 0; i < noc - 1; ++i) {
         ffInsert(row1,i,FF_ONE);
      }
      result |= TestRowIo2(row0,row1,buf,noc);
      free(row0);
      free(row1);
      free(buf);
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
static int TestHdr2(int nor, int noc, PTR buf1, PTR buf2)
{
   FILE *f;
   const char *file_name = "check.1";
   int fld2 = -1, nor2 = -1, noc2 = -1;

   /* Write <buf1> into file
      ---------------------- */
   f = ffWriteHeader(file_name,ffOrder,nor,noc);
   ffWriteRows(f,buf1,nor,noc);
   fclose(f);

   /* Read the file header and check the value.
      ----------------------------------------- */
   memset(buf2,0,ffRowSize(noc) * nor);
   f = ffReadHeader(file_name,&fld2,&nor2,&noc2);
   ASSERT_EQ_INT(fld2, ffOrder);
   ASSERT_EQ_INT(nor2, nor);
   ASSERT_EQ_INT(noc2, noc);

   // Read the rows. 
   ffReadRows(f,buf2,nor,noc);
   fclose(f);

   // Compare <buf1> and <buf2>
   //tstPrintRows("buf1", buf1, nor, noc);
   //tstPrintRows("buf2", buf2, nor, noc);
   ASSERT_EQ_INT(CmpMat(buf1,buf2,nor,noc),0);

   remove(file_name);
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
      long x = 0;

      buf1 = ffAlloc(nor, noc);
      buf2 = ffAlloc(nor, noc);
      for (i = 0; i < nor; ++i) {
         int k;
         PTR p = (PTR)((char *) buf1 + ffSize(i, noc));
         for (k = 0; k < noc; ++k) {
            ffInsert(p,k,ffFromInt(x >> 10));
            x = x * 69069 + 13;
         }
      }

      result |= TestHdr2(nor,noc,buf1,buf2);

      free(buf1);
      free(buf2);
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestSeek2(int nor, int noc, PTR buf1, PTR buf2)
{
   FILE *f;
   const char *file_name = "check.1";
   int fld2, nor2, noc2;
   int i;

   /* Write <buf1> into file
      ---------------------- */
   f = ffWriteHeader(file_name,ffOrder,nor,noc);
   ffWriteRows(f,buf1,nor, noc);
   fclose(f);

   /* Read the rows in reverse order.
      ------------------------------- */
   memset(buf2,0,ffRowSize(noc) * nor);
   f = ffReadHeader(file_name, &fld2, &nor2, &noc2);
   for (i = nor - 1; i >= 0; --i) {
      PTR x = (PTR)((char *)buf2 + ffSize(i, noc));
      ffSeekRow(f, i, noc2);
      ffReadRows(f,x,1,noc);
   }
   fclose(f);

   /* Compare <buf1> and <buf2>
      ------------------------- */
   ASSERT_EQ_INT(CmpMat(buf1,buf2,nor,noc), 0);

   remove(file_name);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Kernel_Seek(int q)
{
   const int bufsize = 100;
   int result = 0;

   for (int noc = 0; result != 0 && noc < 65; ++noc) {
      PTR buf1, buf2;
      int i;
      long x = 0;

      buf1 = ffAlloc(bufsize, noc);
      buf2 = ffAlloc(bufsize, noc);
      for (i = 0; i < bufsize; ++i) {
         int k;
         PTR p = (PTR)((char *)buf1 + ffSize(i, noc));
         for (k = 0; k < noc; ++k) {
            ffInsert(p,k,ffFromInt(x >> 10));
            x = x * 69069 + 13;
         }
      }

      result |= TestSeek2(bufsize, noc,buf1,buf2);

      free(buf1);
      free(buf2);
   }
   return result;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
