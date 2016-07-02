////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for the command line parser
//
// (C) Copyright 1998-2016 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
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
   ASSERT((t = AppGetTextOption(App,"-b --option2",NULL)) != NULL);
   ASSERT(!strcmp(t,"optarg2"));
   ASSERT(!TstHasError());

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
   ASSERT((App = AppAlloc(NULL,ArgC1,ArgV1)) != NULL);
   ASSERT(AppGetArguments(App,2,2) == 2 && !TstHasError());
   ASSERT(AppGetArguments(App,1,1) == -1 && TstHasError());	// too many arguments
   ASSERT(AppGetArguments(App,3,3) == -1 && TstHasError());	// too few arguments
   AppFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

test_F DetectUnknownOption()
{
   TstStartErrorChecking();
   static const char *ArgV1[] = { "---", "-a", "--option1", "--option2" };
   const int ArgC1 = (sizeof(ArgV1) / sizeof(ArgV1[0]));
   MtxApplication_t *App;
   ASSERT((App = AppAlloc(NULL,ArgC1,ArgV1)) != NULL);
   AppGetOption(App,"-a");
   AppGetOption(App,"--option1");
   ASSERT(AppGetArguments(App,0,100) == -1 && TstHasError());
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
   ASSERT((App = AppAlloc(NULL,ArgC1,ArgV1)) != NULL);
   ASSERT(AppGetOption(App,"-a"));
   ASSERT(!AppGetOption(App,"-b") && !TstHasError());	// -b is not an option

   ASSERT(AppGetArguments(App,1,1) == 1 && !TstHasError());
   ASSERT(!strcmp(App->ArgV[0],"-b"));
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
   ASSERT((App = AppAlloc(NULL,argc,argv)) != NULL);
   ASSERT(AppGetIntOption(App,"-a",42,1,10) == 10);
   ASSERT(AppGetIntOption(App,"-b --bbb",42,-20,-19) == -20);
   ASSERT(AppGetIntOption(App,"-c",42,0,-1) == 3);
   ASSERT_EQ_INT(AppGetArguments(App,0,0), 0);
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
   ASSERT((App = AppAlloc(NULL,argc,argv)) != NULL);
   ASSERT(AppGetIntOption(App,"-a",42,1,0) == 42 && TstHasError());	// malformed value
   ASSERT(AppGetIntOption(App,"-b --bbb",42,21,999) == 42 && TstHasError()); // out of range
   ASSERT(AppGetIntOption(App,"-c",42,0,29) == 42 && TstHasError());	// out of range
   ASSERT(AppGetArguments(App,0,0) == 0);
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
   ASSERT((App = AppAlloc(NULL,argc,argv)) != NULL);
   ASSERT(AppGetOption(App,"-a") && !TstHasError());
   ASSERT(AppGetOption(App,"-b") && !TstHasError());
   ASSERT(AppGetArguments(App,0,110) != 0 && TstHasError());
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
   ASSERT((App = AppAlloc(NULL,argc,argv)) != NULL);
   ASSERT(AppGetCountedOption(App,"-a --all") == 4 && !TstHasError());
   ASSERT(AppGetCountedOption(App,"-b") == 3 && !TstHasError());
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
   ASSERT((App = AppAlloc(NULL,argc,argv)) != NULL);
   ASSERT(MtxMessageLevel <= -1000);			// --quiet
   ASSERT(!strcmp(MtxBinDir,"binbinbin"));		// -B
   ASSERT(!strcmp(MtxLibDir,"libliblib"));		// -L
   ASSERT(AppGetArguments(App,0,100) == 0 && !TstHasError());
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
   ASSERT((App = AppAlloc(NULL,argc,argv)) != NULL);
   ASSERT_EQ_INT(MtxMessageLevel, 3);
   MtxMessageLevel = 0;
   ASSERT(!strcmp(MtxBinDir,"BINBIN"));
   ASSERT(!strcmp(MtxLibDir,"LIBLIB"));
   ASSERT(AppGetArguments(App,0,100) == 0 && !TstHasError());
   AppFree(App);
   strcpy(MtxBinDir,".");
   strcpy(MtxLibDir,".");
}

