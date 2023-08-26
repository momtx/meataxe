////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check stfXXX()
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult StructuredTextFile1()
{
   const char *file_name = "check.tmp";
   const char *string1 = "\t this is a\r\tst\a\bri\"ng\f\n   ";
   const int num1 = 42;
   const int vec1[10] = { -1,0,1,2,3,4,5,6,7,8 };
   char string2[100];
   int num2;
   int vec2[10] = { 0 };
   int vec2size = 10;
   StfData *f;

   ASSERT((f = stfOpen(file_name,"w")) != NULL);
   ASSERT(stfWriteValue(f,"StfTest","rec()") == 0);
   ASSERT(stfWriteString(f,"StfTest.String1",string1) == 0);

   ASSERT(stfWriteInt(f,"StfTest.Integer1",num1) == 0);
   ASSERT(stfWriteVector(f,"StfTest.Vector1",10,vec1) == 0);
   stfClose(f);

   ASSERT((f = stfOpen(file_name,"rb")) != NULL);
   ASSERT(stfReadLine(f) == 0);
   ASSERT(strcmp(stfGetName(f),"StfTest") == 0);
      
   ASSERT(stfReadLine(f) == 0);
   ASSERT(strcmp(stfGetName(f),"StfTest.String1") == 0);
   ASSERT(stfGetString(f,string2,sizeof(string2)) == 0);
   ASSERT(strcmp(string1,string2) == 0);
   
   ASSERT(stfReadLine(f) == 0);
   ASSERT(strcmp(stfGetName(f),"StfTest.Integer1") == 0);
   ASSERT(stfGetInt(f,&num2) == 0);
   ASSERT(num1 == num2);
   
   ASSERT(!stfReadLine(f));
   ASSERT(strcmp(stfGetName(f),"StfTest.Vector1") == 0);
   ASSERT(stfGetVector(f,&vec2size,vec2) == 0);
   ASSERT(memcmp(vec1,vec2,sizeof(vec1)) == 0);
   stfClose(f);

   remove(file_name);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult StructuredTextFile2()
{
   const char *file_name = "check.tmp";
   int vec1[1000];
   int vec2[1000];
   int vec2size = 1000;
   int i;
   StfData *f;

   for (i = 0; i < 1000; ++i) {
      vec1[i] = i;
   }
   ASSERT((f = stfOpen(file_name,"w")) != NULL);
   ASSERT(stfWriteVector(f,"StfTest.Vector1",1000,vec1) == 0);
   stfClose(f);

   ASSERT((f = stfOpen(file_name,"rb")) != NULL);
   ASSERT(stfReadLine(f) == 0);
   ASSERT(strcmp(stfGetName(f),"StfTest.Vector1") == 0);
   ASSERT(stfGetVector(f,&vec2size,vec2) == 0);
   ASSERT(memcmp(vec1,vec2,sizeof(vec1)) == 0);
   stfClose(f);

   remove(file_name);
   return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
