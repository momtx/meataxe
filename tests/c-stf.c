////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Check StfXXX()
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F StructuredTextFile1()
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

   if ((f = StfOpen(file_name,FM_CREATE)) == NULL) {
      Error("StfOpen() filed");
   }
   if (StfWriteValue(f,"StfTest","rec()") != 0) {
      Error("StfWriteValue() failed");
   }
   if (StfWriteString(f,"StfTest.String1",string1) != 0) {
      Error("StfWriteString() failed");
   }
   if (StfWriteInt(f,"StfTest.Integer1",num1) != 0) {
      Error("StfWriteInt() failed");
   }
   if (StfWriteVector(f,"StfTest.Vector1",10,vec1) != 0) {
      Error("StfWriteVector() failed");
   }
   StfClose(f);

   if ((f = StfOpen(file_name,FM_READ)) == NULL) {
      Error("StfOpen() filed");
   }
   if (StfReadLine(f) || strcmp(StfGetName(f),"StfTest")) {
      Error("Header not found");
   }
   if (StfReadLine(f) || strcmp(StfGetName(f),"StfTest.String1")
       || StfGetString(f,string2,sizeof(string2)) || strcmp(string1,string2)) {
      Error("Read string failed");
   }
   if (StfReadLine(f) || strcmp(StfGetName(f),"StfTest.Integer1")
       || StfGetInt(f,&num2) || (num1 != num2)) {
      Error("Read integer failed");
   }
   if (StfReadLine(f) || strcmp(StfGetName(f),"StfTest.Vector1")
       || StfGetVector(f,&vec2size,vec2) || memcmp(vec1,vec2,sizeof(vec1))) {
      Error("Read vector failed");
   }
   StfClose(f);

   remove(file_name);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F StructuredTextFile2()
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
   if ((f = StfOpen(file_name,FM_CREATE)) == NULL) {
      Error("StfOpen() filed");
   }
   if (StfWriteVector(f,"StfTest.Vector1",1000,vec1) != 0) {
      Error("StfWriteVector() failed");
   }
   StfClose(f);

   if ((f = StfOpen(file_name,FM_READ)) == NULL) {
      Error("StfOpen() filed");
   }
   if (StfReadLine(f) || strcmp(StfGetName(f),"StfTest.Vector1")
       || StfGetVector(f,&vec2size,vec2) || memcmp(vec1,vec2,sizeof(vec1))) {
      Error("Read vector failed");
   }
   StfClose(f);

   remove(file_name);
}
