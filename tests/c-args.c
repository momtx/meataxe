////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for the argument parser
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "check.h"

#include <string.h>

static int ErrorFlag = 0;
static MtxApplication_t *App;

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Command line parsing

test_F CommandLineParsing1()
{
   TstStartErrorChecking();

   const char *t;
   static const char *ArgV1[] =
   { "---", "-a", "--option1", "--option2", "optarg2", "arg1", "arg2" };
   int ArgC1 = (sizeof(ArgV1) / sizeof(ArgV1[0]));

   if ((App = AppAlloc(NULL,ArgC1,ArgV1)) == NULL) {
      TST_FAIL("AppAlloc() failed");
   }

   if (AppGetOption(App,"--option11")) {
      TST_FAIL("Option --option11 recognized");
   }
   if (!AppGetOption(App,"--option1") || TstHasError()) {
      TST_FAIL("Option --option1 not recognized");
   }
   if (!AppGetOption(App,"-a --aaaa") || TstHasError()) {
      TST_FAIL("Option -a not recognized");
   }
   if (AppGetOption(App,"-a --aaaa") || TstHasError()) {
      TST_FAIL("Option -a repeated");
   }
   t = AppGetTextOption(App,"-b --option2",NULL);
   if ((t == NULL) || strcmp(t,"optarg2") || TstHasError()) {
      TST_FAIL("Text option not recognized");
   }

   AppFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F ArgumentCountChecking()
{
   TstStartErrorChecking();

   static const char *ArgV1[] =
   { "---", "arg1", "arg2" };
   int ArgC1 = (sizeof(ArgV1) / sizeof(ArgV1[0]));

   if ((App = AppAlloc(NULL,ArgC1,ArgV1)) == NULL) {
      TST_FAIL("AppAlloc() failed");
   }

   if ((AppGetArguments(App,2,2) != 2) || TstHasError()) {
      TST_FAIL("AppGetArguments() failed");
   }
   if ((AppGetArguments(App,1,1) != -1) || !TstHasError()) {
      TST_FAIL("AppGetArguments(App,1,1) failed");
   }
   if ((AppGetArguments(App,3,3) != -1) || !TstHasError()) {
      TST_FAIL("AppGetArguments(3,3) failed");
   }
   AppFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F DetectUnknownOption()
{
   TstStartErrorChecking();

   static const char *ArgV1[] =
   { "---", "-a", "--option1", "--option2" };
   int ArgC1 = (sizeof(ArgV1) / sizeof(ArgV1[0]));
   if ((App = AppAlloc(NULL,ArgC1,ArgV1)) == NULL) {
      TST_FAIL("AppAlloc() failed");
   }
   AppGetOption(App,"-a");
   AppGetOption(App,"--option1");
   if ((AppGetArguments(App,0,100) != -1) || !TstHasError()) {
      TST_FAIL("Unknown option not detected");
   }
   AppFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void MyErrorHandler(const MtxErrorRecord_t *err)
{
   ErrorFlag = (err != NULL) ? 1 : 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int CheckError()
{
   int i = ErrorFlag;
   ErrorFlag = 0;
   return i;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/* Test 2
   Erkennung von '--'
 */

static void Test2()
{
   static const char *ArgV1[] =
   { "---", "-a", "--", "-b" };
   int ArgC1 = (sizeof(ArgV1) / sizeof(ArgV1[0]));

   if ((App = AppAlloc(NULL,ArgC1,ArgV1)) == NULL) {
      TST_FAIL("AppAlloc() failed");
   }

   if (!AppGetOption(App,"-a")) {
      TST_FAIL("Option -a not recognized");
   }
   if (AppGetOption(App,"-b") || CheckError()) {
      TST_FAIL("Option -b found after '--'");
   }
   if ((AppGetArguments(App,1,1) != 1) || CheckError()) {
      TST_FAIL("AppGetArguments() failed");
   }
   if (strcmp(App->ArgV[0],"-b")) {
      TST_FAIL("Argument after '--' not found");
   }
   AppFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/* Test 3
   Integer-Optionen
 */

static void Test3()
{
   static const char *argv[] =
   { "-", "-a", "10", "--bbb", "-20", "-c", "3" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   if ((App = AppAlloc(NULL,argc,argv)) == NULL) {
      TST_FAIL("AppAlloc() failed");
   }
   if (AppGetIntOption(App,"-a",42,1,10) != 10) {
      TST_FAIL("Option -a 10 not recognized");
   }
   if (AppGetIntOption(App,"-b --bbb",42,-20,-19) != -20) {
      TST_FAIL("Option -bbb 20 not recognized");
   }
   if (AppGetIntOption(App,"-c",42,0,-1) != 3) {
      TST_FAIL("Option -c 3 not recognized");
   }
   if (AppGetArguments(App,0,0) != 0) {
      TST_FAIL("AppGetArguments() failed");
   }
   AppFree(App);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/* Test 4: Integer option errors */

static void Test4()
{
   static const char *argv[] =
   { "-", "-a", "1x0", "--bbb", "20", "-c", "30" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   if ((App = AppAlloc(NULL,argc,argv)) == NULL) {
      TST_FAIL("AppAlloc() failed");
   }
   if ((AppGetIntOption(App,"-a",42,1,0) != 42) || !CheckError()) {
      TST_FAIL("Error in option '-a 1x0' not found");
   }
   if ((AppGetIntOption(App,"-b --bbb",42,21,999) != 42) || !CheckError()) {
      TST_FAIL("Range check failed");
   }
   if ((AppGetIntOption(App,"-c",42,0,29) != 42) || !CheckError()) {
      TST_FAIL("Range check 2 failed");
   }
   if (AppGetArguments(App,0,0) != 0) {
      TST_FAIL("AppGetArguments() failed");
   }
   AppFree(App);
}


/* Test 5: Option after Argument */

static void Test5()
{
   static const char *argv[] =
   { "-", "-a", "xxx", "-b", "yyy" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   if ((App = AppAlloc(NULL,argc,argv)) == NULL) {
      TST_FAIL("AppAlloc() failed");
   }
   if (!AppGetOption(App,"-a") || CheckError()) {
      TST_FAIL("Option '-a' not recognized");
   }
   if (!AppGetOption(App,"-b") || CheckError()) {
      TST_FAIL("Option '-b' not recognized");
   }
   if ((AppGetArguments(App,0,110) == 0) || !CheckError()) {
      TST_FAIL("AppGetArguments())=0 on invalid data");
   }
   AppFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/* Test 6: Counted Option */

static void Test6()
{
   static const char *argv[] =
   { "-", "-a", "-b", "-a", "-bb", "--all", "-ca" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   if ((App = AppAlloc(NULL,argc,argv)) == NULL) {
      TST_FAIL("AppAlloc() failed");
   }
   if ((AppGetCountedOption(App,"-a --all") != 4) || CheckError()) {
      TST_FAIL("Option '--all' not recognized");
   }
   if ((AppGetCountedOption(App,"-b") != 3) || CheckError()) {
      TST_FAIL("Option '-b' not recognized");
   }
   AppFree(App);
}


/* Test 7: Common options */

static void Test7()
{
   static const char *argv[] =
   { "-", "--quiet", "-B", "binbinbin", "-L", "libliblib" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   if ((App = AppAlloc(NULL,argc,argv)) == NULL) {
      TST_FAIL("AppAlloc() failed");
   }
   if (MtxMessageLevel > -1000) {
      TST_FAIL("Option --quiet not processed");
   }
   if (strcmp(MtxBinDir,"binbinbin")) {
      TST_FAIL("Option -B not processed");
   }
   if (strcmp(MtxLibDir,"libliblib")) {
      TST_FAIL("Option -L not processed");
   }
   if ((AppGetArguments(App,0,100) != 0) || CheckError()) {
      TST_FAIL("Common options missed");
   }
   AppFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/* Test 8: Common options (b) */

static void Test8()
{
   static const char *argv[] =
   { "-", "-V", "-VV", "--mtxlib", "LIBLIB", "--mtxbin", "BINBIN" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   if ((App = AppAlloc(NULL,argc,argv)) == NULL) {
      TST_FAIL("AppAlloc() failed");
   }
   if (MtxMessageLevel != 3) {
      TST_FAIL("Option -V not processed correctly");
   }
   MtxMessageLevel = 0;
   if (strcmp(MtxBinDir,"BINBIN")) {
      TST_FAIL("Option --mtxbin not processed");
   }
   if (strcmp(MtxLibDir,"LIBLIB")) {
      TST_FAIL("Option --mtxlib not processed");
   }
   if ((AppGetArguments(App,0,100) != 0) || CheckError()) {
      TST_FAIL("Common options missed");
   }
   AppFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F CommandLineProcessing()
{
   MtxErrorHandler_t *oldhandler;
   oldhandler = MtxSetErrorHandler(MyErrorHandler);
   Test2();
   Test3();
   Test4();
   Test5();
   Test6();
   Test7();
   Test8();

   MtxSetErrorHandler(oldhandler);
   strcpy(MtxBinDir,".");
   strcpy(MtxLibDir,".");
}
