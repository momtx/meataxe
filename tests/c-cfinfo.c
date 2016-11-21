////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for .cfinfo related functions
//
// (C) Copyright 2016 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F CfinfoFileNotFound()
{
   TstStartErrorChecking();

   Lat_Info info;
   ASSERT(Lat_ReadInfo(&info, "/file_not_found") == -1);
   ASSERT(TstHasError());
}

////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct Tst_TempFile {
   char *base_name;
   char *name;
   struct Tst_TempFile* next;
} Tst_TempFile;
static Tst_TempFile *TempFiles = NULL;
static int TempFileId = 0;

const Tst_TempFile* TstCreateTemporaryFile(const char *ext, const char *data)
{
   ++TempFileId;
   Tst_TempFile *tf = (Tst_TempFile*) malloc(sizeof(Tst_TempFile));
   tf->next = TempFiles;
   TempFiles = tf;

   // base name
   const size_t bfnl = 20;
   tf->base_name = (char*) malloc(bfnl);
   snprintf(tf->base_name,bfnl,"tmp%08x",TempFileId);

   // full name
   const size_t fnl = bfnl + strlen(ext);
   tf->name = (char*) malloc(fnl);
   snprintf(tf->name,fnl,"%s%s",tf->base_name,ext);

   // create file
   const int fd = open(tf->name,O_WRONLY|O_CREAT|O_TRUNC,0600);
   write(fd,data,strlen(data));
   close(fd);

   return tf;
}

void Tst_RemoveTempFiles()
{
   while (TempFiles) {
      Tst_TempFile *tf = TempFiles;
      TempFiles = TempFiles->next;
      remove(tf->name);
      free(tf->name);
      free(tf->base_name);
      free(tf);
   }
}

test_F CfinfoBadIdWord()
{
   static const char FILE_DATA[] =
      "CFInfo := rec();\n"
      "CFInfo.NGen := 2;\n"
      "CFInfo.Field := 2;\n"
      "CFInfo.NCF := 1;\n"
      "CFInfo.ConstituentNames := [\"x\"];\n"
      "CFInfo.Dimension := [10];\n"
      "CFInfo.Number := [0];\n"
      "CFInfo.Multiplicity := [1];\n"
      "CFInfo.SplittingField := [1];\n"
      "CFInfo.NMountains := [0];\n"
      "CFInfo.NDottedLines := [0];\n"
      "CFInfo.PeakWord := [17,2,1,1,1];\n"      /// syntax error
      "CFInfo.IdWord := [[3,2,1,0,1]];\n"
      "CFInfo.NSocles := 0;\n"
      "CFInfo.Socles := [];\n"
      "CFInfo.NHeads := 0;\n"
      "CFInfo.Heads := [];\n";

   TstStartErrorChecking();
   const Tst_TempFile* const tf = TstCreateTemporaryFile(".cfinfo",FILE_DATA);

   Lat_Info info;
   ASSERT(Lat_ReadInfo(&info, tf->base_name) == -1);
   ASSERT(TstHasError());
   ASSERT(info.NCf == 1);
   //Tst_RemoveTempFiles();
}


#if 0

   const char *file_name = "check.tmp";
   const char *string1 = "\t this is a\r\tst\a\bri\"ng\f\n   ";
   const int num1 = 42;
   const int vec1[10] = { -1,0,1,2,3,4,5,6,7,8 };
   char string2[100];
   int num2;
   int vec2[10] = { 0 };
   int vec2size = 10;
   StfData *f;

   if ((f = StfOpen(file_name,FM_CREATE)) == NULL) {
      TST_FAIL("StfOpen() filed");
   }
   if (StfWriteValue(f,"StfTest","rec()") != 0) {
      TST_FAIL("StfWriteValue() failed");
   }
   if (StfWriteString(f,"StfTest.String1",string1) != 0) {
      TST_FAIL("StfWriteString() failed");
   }
   if (StfWriteInt(f,"StfTest.Integer1",num1) != 0) {
      TST_FAIL("StfWriteInt() failed");
   }
   if (StfWriteVector(f,"StfTest.Vector1",10,vec1) != 0) {
      TST_FAIL("StfWriteVector() failed");
   }
   StfClose(f);

   if ((f = StfOpen(file_name,FM_READ)) == NULL) {
      TST_FAIL("StfOpen() filed");
   }
   if (StfReadLine(f) || strcmp(StfGetName(f),"StfTest")) {
      TST_FAIL("Header not found");
   }
   if (StfReadLine(f) || strcmp(StfGetName(f),"StfTest.String1")
       || StfGetString(f,string2,sizeof(string2)) || strcmp(string1,string2)) {
      TST_FAIL("Read string failed");
   }
   if (StfReadLine(f) || strcmp(StfGetName(f),"StfTest.Integer1")
       || StfGetInt(f,&num2) || (num1 != num2)) {
      TST_FAIL("Read integer failed");
   }
   if (StfReadLine(f) || strcmp(StfGetName(f),"StfTest.Vector1")
       || StfGetVector(f,&vec2size,vec2) || memcmp(vec1,vec2,sizeof(vec1))) {
      TST_FAIL("Read vector failed");
   }
   StfClose(f);

   remove(file_name);
}
#endif

// vim:sw=3:et
