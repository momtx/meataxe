////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for the command line parser
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>
#include "testing.h"

#include <string.h>


////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_LongOptionArgumentCannotBeSeparated()
{
   static char* const argv[] = { "---", "--option", "optval" };
   const int argc = (sizeof(argv) / sizeof(argv[0]));
   MtxApplication_t* app = appAlloc(NULL, argc, argv);

   // "optval is not recognized as argument to "--option".
   ASSERT_EQ_STRING(appGetTextOption(app,"-o --option", "dflt"), "dflt");
   // "optval is treated as normal argument.
   ASSERT_EQ_INT(appGetArguments(app, 0, 1), 1);
   ASSERT_EQ_STRING(app->argV[0], "optval");

   appFree(app);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_FailsOnMissingMandatoryArgumentOfLongOption()
{
   static char* const argv[] =
   { "---", "--option", "optval" };
   const int argc = (sizeof(argv) / sizeof(argv[0]));
   MtxApplication_t* app = appAlloc(NULL, argc, argv);

   ASSERT_ABORT(appGetTextOption(app,"-o --option", NULL));

   appFree(app);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_RecognizedLongOptionIsConsumed()
{
   static char* const argv[] =
   { "---", "--option1", "-ABC" };
   const int argc = (sizeof(argv) / sizeof(argv[0]));
   MtxApplication_t* app = appAlloc(NULL, argc, argv);

   // Read all options
   ASSERT_EQ_INT(appGetOption(app,"-A"), 1);
   ASSERT_EQ_INT(appGetOption(app,"-o --option1"), 1);
   ASSERT_EQ_INT(appGetOption(app,"-C"), 1);
   ASSERT_EQ_INT(appGetOption(app,"-B"), 1);

   // Try again -> not found
   ASSERT_EQ_INT(appGetOption(app,"-A"), 0);
   ASSERT_EQ_INT(appGetOption(app,"-o --option1"), 0);
   ASSERT_EQ_INT(appGetOption(app,"-C"), 0);
   ASSERT_EQ_INT(appGetOption(app,"-B"), 0);

   appFree(app);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_LongOptionWithValue()
{
   static char* argv[] = { "---", "--option1=value1", "--option2=value2", "arg" };
   int argc = (sizeof(argv) / sizeof(argv[0]));
   MtxApplication_t* app = appAlloc(NULL, argc, argv);

   ASSERT_EQ_STRING(appGetTextOption(app,"--option1", NULL), "value1");
   ASSERT_EQ_STRING(appGetTextOption(app,"--option2", "dflt2"), "value2");
   ASSERT_EQ_INT(appGetArguments(app, 1, 1), 1);
   appFree(app);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_LongOptionWithoutValueWithDefault()
{
   static char* argv[] = { "---", "--option1", "arg" };
   static const int argc = (sizeof(argv) / sizeof(argv[0]));
   MtxApplication_t* app = appAlloc(NULL, argc, argv);
   ASSERT_EQ_STRING(appGetTextOption(app,"--option1", "dflt"),"dflt");
   appFree(app);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_LongOptionWithoutValueWithoutDefault()
{
   static char* argv[] = { "---", "--option1", "arg" };
   static const int argc = (sizeof(argv) / sizeof(argv[0]));
   MtxApplication_t* app = appAlloc(NULL, argc, argv);
   ASSERT_ABORT(appGetTextOption(app,"--option1", NULL));
   appFree(app);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_LongOptionNotPresent()
{
   static char* argv[] = { "---", "arg" };
   static const int argc = (sizeof(argv) / sizeof(argv[0]));
   MtxApplication_t* app = appAlloc(NULL, argc, argv);
   ASSERT_EQ_STRING(appGetTextOption(app,"--option1", NULL), NULL);
   ASSERT_EQ_STRING(appGetTextOption(app,"--option2", "dflt"), NULL);
   ASSERT_EQ_INT(appGetArguments(app, 1, 1), 1);
   appFree(app);
   return 0;
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
   ASSERT(!strcmp(App->argV[0],"-b"));
   appFree(App);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_IntegerOptions()
{
   static char *argv[] =
   { "-", "-a", "10", "--bbb=-20", "-c", "3" };
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

TstResult App_ShortOptionsWithValueCannotBeMerged()
{
   static char *argv[] = { "-", "-ab", "aval" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   MtxApplication_t *App;
   ASSERT((App = appAlloc(NULL,argc,argv)) != NULL);
   ASSERT_EQ_INT(appGetOption(App,"-b"), 1);
   ASSERT_ABORT(appGetIntOption(App,"-a",0,0,100));
   appFree(App);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_IntegerOptionErrorHandling()
{
   static char *argv[] =
   { "-", "-a", "1x0", "--bbb=20", "-c", "30" };
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

//TODO: rewrite using appGetOption
//TstResult App_CountedOption()
//{
//   static char *argv[] =
//   { "-", "-a", "-b", "-a", "-bb", "--all", "-ca" };
//   int argc = (sizeof(argv) / sizeof(argv[0]));
//
//   MtxApplication_t *App;
//   ASSERT((App = appAlloc(NULL,argc,argv)) != NULL);
//   ASSERT(appGetCountedOption(App,"-a --all") == 4);
//   ASSERT(appGetCountedOption(App,"-b") == 3);
//   appFree(App);
//   return 0;
//}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_CommonOptions()
{
   static char *argv[] =
   { "-", "--quiet", "-L", "libliblib" };
   int argc = (sizeof(argv) / sizeof(argv[0]));

   mtxCleanupLibrary();
   MtxApplication_t *App;
   ASSERT((App = appAlloc(NULL,argc,argv)) != NULL);
   ASSERT_EQ_INT(logGetDefaultThreshold(), MTX_LOG_WARNING);   // --quiet
   ASSERT(!strcmp(mtxLibraryDirectory(),"libliblib"));	// -L
   ASSERT(appGetArguments(App,0,100) == 0);
   appFree(App);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult App_CommonOptions2()
{
   static char *argv[] = { "-", "-V", "--verbose", "--mtxlib=LIBLIB" };
   int argc = (sizeof(argv) / sizeof(argv[0]));
   mtxCleanupLibrary();
   MtxApplication_t *App;

   ASSERT((App = appAlloc(NULL,argc,argv)) != NULL);

   ASSERT_EQ_INT(logGetDefaultThreshold(), MTX_LOG_DEBUG2);
   ASSERT(!strcmp(mtxLibraryDirectory(),"LIBLIB"));
   ASSERT(appGetArguments(App,0,100) == 0);
   appFree(App);
   return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
