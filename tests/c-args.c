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


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Command line parsing

test_F CommandLineParsing1()
{
   TstStartErrorChecking();

   const char *t;
   static const char *ArgV1[] =
   { "---", "-a", "--option1", "--option2", "optarg2", "arg1", "arg2" };
   int ArgC1 = (sizeof(ArgV1) / sizeof(ArgV1[0]));

   MtxApplication_t *App;
   ASSERT((App = AppAlloc(NULL,ArgC1,ArgV1)) != NULL);

   ASSERT(!AppGetOption(App,"--option11") && !TstHasError());	// not present
   ASSERT(AppGetOption(App,"--option1") && !TstHasError());
   ASSERT(AppGetOption(App,"-a --aaaa") && !TstHasError());	// short and long name
   ASSERT(!AppGetOption(App,"-a --aaaa") && !TstHasError());	// repeated
  
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

   MtxApplication_t *App;
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
   MtxApplication_t *App;
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

/// Using "--" to designate the end of options

test_F DoubleDash()
{
   TstStartErrorChecking();

   static const char *ArgV1[] =
   { "---", "-a", "--", "-b" };
   int ArgC1 = (sizeof(ArgV1) / sizeof(ArgV1[0]));

   MtxApplication_t *App;
   if ((App = AppAlloc(NULL,ArgC1,ArgV1)) == NULL) {
      TST_FAIL("AppAlloc() failed");
   }

   if (!AppGetOption(App,"-a")) {
      TST_FAIL("Option -a not recognized");
   }
   if (AppGetOption(App,"-b") || TstHasError()) {
      TST_FAIL("Option -b found after '--'");
   }
   if ((AppGetArguments(App,1,1) != 1) || TstHasError()) {
      TST_FAIL("AppGetArguments() failed");
   }
   if (strcmp(App->ArgV[0],"-b")) {
      TST_FAIL("Argument after '--' not found");
   }
   AppFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F IntegerOptions()
{
   TstStartErrorChecking();

   static const char *argv[] =
   { "-", "-a", "10", "--bbb", "-20", "-c", "3" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   MtxApplication_t *App;
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

test_F IntegerOptionErrorHandling()
{
   TstStartErrorChecking();
   static const char *argv[] =
   { "-", "-a", "1x0", "--bbb", "20", "-c", "30" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   MtxApplication_t *App;
   if ((App = AppAlloc(NULL,argc,argv)) == NULL) {
      TST_FAIL("AppAlloc() failed");
   }
   if ((AppGetIntOption(App,"-a",42,1,0) != 42) || !TstHasError()) {
      TST_FAIL("Error in option '-a 1x0' not found");
   }
   if ((AppGetIntOption(App,"-b --bbb",42,21,999) != 42) || !TstHasError()) {
      TST_FAIL("Range check failed");
   }
   if ((AppGetIntOption(App,"-c",42,0,29) != 42) || !TstHasError()) {
      TST_FAIL("Range check 2 failed");
   }
   if (AppGetArguments(App,0,0) != 0) {
      TST_FAIL("AppGetArguments() failed");
   }
   AppFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F OptionAfterArgument()
{
   TstStartErrorChecking();
   static const char *argv[] =
   { "-", "-a", "xxx", "-b", "yyy" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   MtxApplication_t *App;
   if ((App = AppAlloc(NULL,argc,argv)) == NULL) {
      TST_FAIL("AppAlloc() failed");
   }
   if (!AppGetOption(App,"-a") || TstHasError()) {
      TST_FAIL("Option '-a' not recognized");
   }
   if (!AppGetOption(App,"-b") || TstHasError()) {
      TST_FAIL("Option '-b' not recognized");
   }
   if ((AppGetArguments(App,0,110) == 0) || !TstHasError()) {
      TST_FAIL("AppGetArguments())=0 on invalid data");
   }
   AppFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F CountedOption()
{
   TstStartErrorChecking();
   static const char *argv[] =
   { "-", "-a", "-b", "-a", "-bb", "--all", "-ca" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   MtxApplication_t *App;
   if ((App = AppAlloc(NULL,argc,argv)) == NULL) {
      TST_FAIL("AppAlloc() failed");
   }
   if ((AppGetCountedOption(App,"-a --all") != 4) || TstHasError()) {
      TST_FAIL("Option '--all' not recognized");
   }
   if ((AppGetCountedOption(App,"-b") != 3) || TstHasError()) {
      TST_FAIL("Option '-b' not recognized");
   }
   AppFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F CommonOptions()
{
   TstStartErrorChecking();
   static const char *argv[] =
   { "-", "--quiet", "-B", "binbinbin", "-L", "libliblib" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   MtxApplication_t *App;
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
   if ((AppGetArguments(App,0,100) != 0) || TstHasError()) {
      TST_FAIL("Common options missed");
   }
   AppFree(App);
   strcpy(MtxBinDir,".");
   strcpy(MtxLibDir,".");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F CommonOptions2()
{
   TstStartErrorChecking();

   static const char *argv[] =
   { "-", "-V", "-VV", "--mtxlib", "LIBLIB", "--mtxbin", "BINBIN" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   MtxApplication_t *App;
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
   if ((AppGetArguments(App,0,100) != 0) || TstHasError()) {
      TST_FAIL("Common options missed");
   }
   AppFree(App);
   strcpy(MtxBinDir,".");
   strcpy(MtxLibDir,".");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F CommandLineProcessing()
{
}
