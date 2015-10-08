////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check functions for the OS interface.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"

#include <string.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static void CheckMem(char *x, char val, int len)
{
   while (len > 0 && *x == val) {
      --len;
      ++x;
   }
   if (len > 0) {
      Error("CheckMem(val=%d, len=%d) failed",(int)val,len);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestAlloc()
{
   const int nblk = 10;
   int i;
   char *x[10];

   for (i = 0; i < nblk; ++i) {
      if ((x[i] = SysMalloc(0)) == NULL) {
         Error("SysMalloc(0) failed");
      }
      if ((x[i] = SysRealloc(x[i],100)) == NULL) {
         Error("SysRealloc(100) failed");
      }
   }
   for (i = 0; i < nblk; ++i) {
      memset(x[i],33,100);
   }
   for (i = 0; i < nblk; ++i) {
      CheckMem(x[i],33,100);
   }
   for (i = 0; i < nblk; ++i) {
      if ((x[i] = SysRealloc(x[i],i * 20)) == NULL) {
         Error("SysRealloc(100) failed");
      }
   }
   for (i = 0; i < nblk; ++i) {
      memset(x[i],44,i * 20);
   }
   for (i = 0; i < nblk; ++i) {
      CheckMem(x[i],44,i * 20);
   }
   for (i = 0; i < nblk; ++i) {
      SysFree(x[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

//FILE *SysFopen(const char *name, int mode);
//int SysFseek(FILE *file, long pos);

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestFiles()
{
   FILE *f;
   static char Text[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVEXYZ";
   static char Text1[] = "0123401234ABCDEFGHIJKLMNOPQRSTUVEXYZ";
   static char Text2[] = "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx";

   f = SysFopen("__@@$$xsk",FM_READ | FM_NOERROR);
   if (f != NULL) {
      Error("Non-existent file opened successfully");
   }
   f = SysFopen("check1",FM_CREATE);
   if (f == NULL) {
      Error("Cannot create file");
   }
   fwrite(Text,1,10,f);
   SysFseek(f,5);
   fwrite(Text,1,5,f);
   fclose(f);
   f = SysFopen("check1",FM_APPEND);
   if (f == NULL) {
      Error("Cannot open file");
   }
   fwrite(Text + 10,1,26,f);
   fclose(f);

   f = SysFopen("check1",FM_READ);
   if (f == NULL) {
      Error("Cannot open file");
   }
   fread(Text2,1,36,f);
   fclose(f);

   if (memcmp(Text1,Text2,sizeof(Text)) != 0) {
      Error("Compare error");
   }
   remove("check1");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F OsFunctions()
{
   TestAlloc();
   TestFiles();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/* Create values between -2^31 and 2^31 */
#define VAL(i) ((long) ((69069 * (i) + 1) & 0x7FFFFFFF) * (long) (1 - 2 * (i % 2)))

static void TestIntIo1(long *buf, int bufsize, int safe)
{
   int i, k;
   FILE *f;

   for (i = 0; i < bufsize; ++i) {
      buf[i] = VAL(i);
   }
   f = SysFopen("check1",FM_CREATE);
   for (i = k = 0; k < bufsize; ++i) {
      int n = (i <= bufsize - k) ? i : bufsize - k;
      int nwritten = SysWriteLong(f,buf + k,n);
      if (nwritten != n) {
         Error("SysWriteLong(%d)=%d",n,nwritten);
      }
      k += n;
   }
   fclose(f);

   for (i = 1; i < bufsize; i += i / 10 + 1) {
      f = SysFopen("check1",FM_READ);
      for (k = 0; k < bufsize; ) {
         int to_read = (k + i >= bufsize + safe) ? bufsize + safe - k : i;
         int nr = SysReadLong(f,buf + k,to_read);
         if (nr < 0) {
            Error("SysReadLong() failed");
         }
         if (k + nr > bufsize) {
            Error("Read %d bytes (max %d expected)",k + nr,bufsize);
         }
         k += nr;
      }
      fclose(f);
      for (k = 0; k < bufsize; ++k) {
         if (buf[k] != VAL(k)) {
            Error("Unexpected value %ld at position %d",buf[k],k);
         }
      }
   }

   remove("check1");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void TestIntIo2()
{
   FILE *f;
   char buf[16] = {1,0,0,0, 0,2,0,0, 0,0,3,0, 0,0,0,4};
   long rb[4];

   f = SysFopen("check1",FM_CREATE);
   fwrite(buf,1,16,f);
   fclose(f);
   f = SysFopen("check1",FM_READ);
   if (SysReadLong(f,rb,5) != 4) {
      Error("Read error");
   }
   fclose(f);
   if (rb[0] != 0x00000001) {
      Error("Read %lx, expected 0x00000001",buf[0]);
   }
   if (rb[1] != 0x00000200) {
      Error("Read %lx, expected 0x00000200",buf[0]);
   }
   if (rb[2] != 0x00030000) {
      Error("Read %lx, expected 0x00030000",buf[0]);
   }
   if (rb[3] != 0x04000000) {
      Error("Read %lx, expected 0x04000000",buf[0]);
   }
   remove("check1");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F IntegerIo()
{
   long *buf;
   buf = NALLOC(long,10000 + 2000);
   TestIntIo1(buf,10000,2000);
   TestIntIo2();
   free(buf);
}
