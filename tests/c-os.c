////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check functions for the OS interface.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <string.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static int CheckMem(char *x, char val, int len)
{
   while (len > 0 && *x == val) {
      --len;
      ++x;
   }
   if (len > 0) {
      TST_FAIL("CheckMem(val=%d, len=%d) failed",(int)val,len);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult OS_MallocWithLengthZeroIsNotNull()
{
   char *x = sysMalloc(0);
   ASSERT(x != NULL);
   x = sysRealloc(x, 0);
   ASSERT(x != NULL);
   sysFree(x);

   x = sysRealloc(sysMalloc(100), 0);
   ASSERT(x != NULL);
   sysFree(x);
   
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult OS_sysMalloc_InitializesMemoryWithZero()
{
   const int SIZE = 1000000;
   int result = 0;
   for (int i = 0; result == 0 && i < 10; ++i) {
      int *buf = (int*) sysMalloc(SIZE * sizeof(int));
      int k;
      for (k = 0; k < SIZE && buf[k] == 0; ++k);
      if (k < SIZE) {
         tstFail(__FILE__, __LINE__, __func__, "Memory was not initialized");
         result = 1;
      }
      memset(buf, 0xaa, SIZE * sizeof(int));
      sysFree(buf);
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult OS_Alloc()
{
   const int nblk = 10;
   int i;
   char *x[10];
   int result = 0;

   for (i = 0; i < nblk; ++i) {
      ASSERT((x[i] = sysMalloc(100)) != NULL);
   }
   for (i = 0; i < nblk; ++i) {
      memset(x[i],33,100);
   }
   for (i = 0; result == 0 && i < nblk; ++i) {
      result |= CheckMem(x[i],33,100);
   }
   for (i = 0; i < nblk; ++i) {
      ASSERT((x[i] = sysRealloc(x[i],i * 20)) != NULL);
   }
   for (i = 0; result == 0 && i < nblk; ++i) {
      memset(x[i],44,i * 20);
   }
   for (i = 0; i < nblk; ++i) {
      result |= CheckMem(x[i],44,i * 20);
   }
   for (i = 0; i < nblk; ++i) {
      sysFree(x[i]);
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Os_FileIo()
{
   FILE *f;
   static char Text[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVEXYZ";
   static char Text1[] = "0123401234ABCDEFGHIJKLMNOPQRSTUVEXYZ";
   static char Text2[] = "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx";

   ASSERT((f = sysFopen("__@@$$xsk", "rb::noerror")) == NULL);
   ASSERT((f = sysFopen("check1", "wb")) != NULL);
   fwrite(Text,1,10,f);
   sysFseek(f,5);
   fwrite(Text,1,5,f);
   fclose(f);
   ASSERT((f = sysFopen("check1","a")) != NULL);
   fwrite(Text + 10,1,26,f);
   fclose(f);

   ASSERT((f = sysFopen("check1","rb")) != NULL);
   fread(Text2,1,36,f);
   fclose(f);

   ASSERT(memcmp(Text1,Text2,sizeof(Text)) == 0);
   remove("check1");
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/* Create values between -2^31 and 2^31 */
#define VAL(i) ((long) ((69069 * (i) + 1) & 0x7FFFFFFF) * (long) (1 - 2 * (i % 2)))

static int TestIntIo1(uint32_t *buf, int bufsize)
{
   int i, k;
   FILE *f;

   for (i = 0; i < bufsize; ++i) {
      buf[i] = VAL(i);
   }
   f = sysFopen("check1","wb");
   for (i = k = 0; k < bufsize; ++i) {
      int n = (i <= bufsize - k) ? i : bufsize - k;
      sysWrite32(f, buf + k, n);
      k += n;
   }
   fclose(f);

   for (i = 1; i < bufsize; i += i / 10 + 1) {
      f = sysFopen("check1","rb");
      for (k = 0; k < bufsize; ) {
         int to_read = (k + i >= bufsize) ? bufsize - k : i;
         sysRead32(f,buf + k,to_read);
         k += to_read;
      }
      fclose(f);
      for (k = 0; k < bufsize; ++k) {
         ASSERT_EQ_INT(buf[k], VAL(k));
      }
   }

   remove("check1");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int TestIntIo2()
{
   FILE *f;
   char buf[16] = {1,0,0,0, 0,2,0,0, 0,0,3,0, 0,0,0,4};
   uint32_t rb[4];

   f = sysFopen("check1","wb");
   fwrite(buf,1,16,f);
   fclose(f);
   f = sysFopen("check1","rb");
   sysRead32(f,rb,4);
   fclose(f);
   ASSERT(rb[0] == 0x00000001);
   ASSERT(rb[1] == 0x00000200);
   ASSERT(rb[2] == 0x00030000);
   ASSERT(rb[3] == 0x04000000);
   remove("check1");
   return 0;
}

// TODO: ccheck exception on EOF

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Os_IntegerIo()
{
   uint32_t *buf = NALLOC(uint32_t, 10000);
   int result = TestIntIo1(buf,10000);
   result |= TestIntIo2();
   free(buf);
   return result;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
