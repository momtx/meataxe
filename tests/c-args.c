////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for the command line parser
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>
#include "testing.h"

#include <string.h>


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Command line parsing

static int testOptions(MtxApplication_t *app)
{
   ASSERT(!appGetOption(app,"--option11"));	// not present
   ASSERT(appGetOption(app,"--option1"));
   ASSERT(appGetOption(app,"-a --aaaa"));	// short and long name
   ASSERT(!appGetOption(app,"-a --aaaa"));	// repeated

   const char *t;
   ASSERT((t = appGetTextOption(app,"-b --option2",NULL)) != NULL);
   ASSERT(!strcmp(t,"optarg2"));
   return 0;
}

TstResult App_CommandLineParser()
{
    int result = 0;

   static char *ArgV1[] =
   { "---", "-a", "--option1", "--option2", "optarg2", "arg1", "arg2" };
   int ArgC1 = (sizeof(ArgV1) / sizeof(ArgV1[0]));
   MtxApplication_t *App = appAlloc(NULL, ArgC1, ArgV1);
   result |= testOptions(App);
   appFree(App);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_CanCheckArgumentCount()
{
   static char *ArgV1[] =
   { "---", "arg1", "arg2" };
   int ArgC1 = (sizeof(ArgV1) / sizeof(ArgV1[0]));

   MtxApplication_t *App;
   ASSERT((App = appAlloc(NULL,ArgC1, ArgV1)) != NULL);
   ASSERT(appGetArguments(App,2,2) == 2);
   ASSERT_ABORT(appGetArguments(App,1,1));	// too many arguments
   ASSERT_ABORT(appGetArguments(App,3,3));	// too few arguments
   appFree(App);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_DetectUnknownOption()
{
   static char *ArgV1[] = { "---", "-a", "--option1", "--option2" };
   const int ArgC1 = (sizeof(ArgV1) / sizeof(ArgV1[0]));
   MtxApplication_t *App;
   ASSERT((App = appAlloc(NULL,ArgC1,ArgV1)) != NULL);
   appGetOption(App,"-a");
   appGetOption(App,"--option1");
   ASSERT_ABORT(appGetArguments(App,0,100));
   appFree(App);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Using "--" to designate the end of options

TstResult App_DoubleDash()
{
   static char *ArgV1[] =
   { "---", "-a", "--", "-b" };
   int ArgC1 = (sizeof(ArgV1) / sizeof(ArgV1[0]));

   MtxApplication_t *App;
   ASSERT((App = appAlloc(NULL,ArgC1,ArgV1)) != NULL);
   ASSERT(appGetOption(App,"-a"));
   ASSERT(!appGetOption(App,"-b"));	// -b is not an option

   ASSERT(appGetArguments(App,1,1) == 1);
   ASSERT(!strcmp(App->ArgV[0],"-b"));
   appFree(App);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_IntegerOptions()
{
   static char *argv[] =
   { "-", "-a", "10", "--bbb", "-20", "-c", "3" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   MtxApplication_t *App;
   ASSERT((App = appAlloc(NULL,argc,argv)) != NULL);
   ASSERT(appGetIntOption(App,"-a",42,1,10) == 10);
   ASSERT(appGetIntOption(App,"-b --bbb",42,-20,-19) == -20);
   ASSERT(appGetIntOption(App,"-c",42,0,-1) == 3);
   ASSERT_EQ_INT(appGetArguments(App,0,0), 0);
   appFree(App);
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/* Test 4: Integer option errors */

TstResult App_IntegerOptionErrorHandling()
{
   static char *argv[] =
   { "-", "-a", "1x0", "--bbb", "20", "-c", "30" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   MtxApplication_t *App;
   ASSERT((App = appAlloc(NULL,argc,argv)) != NULL);
   ASSERT_ABORT(appGetIntOption(App,"-a",42,1,0));	// malformed value
   ASSERT_ABORT(appGetIntOption(App,"-b --bbb",42,21,999)); // out of range
   ASSERT_ABORT(appGetIntOption(App,"-c",42,0,29));	// out of range
   ASSERT(appGetArguments(App,0,0) == 0);
   appFree(App);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_OptionAfterArgument()
{
   static char *argv[] =
   { "-", "-a", "xxx", "-b", "yyy" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   MtxApplication_t *App;
   ASSERT((App = appAlloc(NULL,argc,argv)) != NULL);
   ASSERT(appGetOption(App,"-a"));
   ASSERT(appGetOption(App,"-b"));
   ASSERT_ABORT(appGetArguments(App,0,110));
   appFree(App);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_CountedOption()
{
   static char *argv[] =
   { "-", "-a", "-b", "-a", "-bb", "--all", "-ca" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   MtxApplication_t *App;
   ASSERT((App = appAlloc(NULL,argc,argv)) != NULL);
   ASSERT(appGetCountedOption(App,"-a --all") == 4);
   ASSERT(appGetCountedOption(App,"-b") == 3);
   appFree(App);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_CommonOptions()
{
   static char *argv[] =
   { "-", "--quiet", "-L", "libliblib" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   MtxApplication_t *App;
   ASSERT((App = appAlloc(NULL,argc,argv)) != NULL);
   ASSERT(MtxMessageLevel == -1);			// --quiet
   ASSERT(!strcmp(MtxLibDir,"libliblib"));		// -L
   ASSERT(appGetArguments(App,0,100) == 0);
   appFree(App);
   strcpy(MtxLibDir,".");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult CommonOptions2()
{
   static char *argv[] =
   { "-", "-V", "-VV", "--mtxlib", "LIBLIB" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   MtxApplication_t *App;
   ASSERT((App = appAlloc(NULL,argc,argv)) != NULL);
   ASSERT_EQ_INT(MtxMessageLevel, 3);
   MtxMessageLevel = 0;
   ASSERT(!strcmp(MtxLibDir,"LIBLIB"));
   ASSERT(appGetArguments(App,0,100) == 0);
   appFree(App);
   strcpy(MtxLibDir,".");
   return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
