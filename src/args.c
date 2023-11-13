////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Command line handling and application utilities
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Global data

#define MTX_MAX_ARGS 150

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Local data

#define IS_DONE(app,i) (app->isDone[i] == 0xFFFFFFFF)
#define IS_DONE_1(app,i,k) ((app->isDone[i] & (1 << (k))) != 0)
#define MARK_DONE(app,i) (app->isDone[i] = 0xFFFFFFFF)
#define MARK_DONE_1(app,i,k) (app->isDone[i] |= (1 << (k)))

/// @defgroup app Application Framework
/// @{
/// The MeatAxe library provides a minimal framework for applications. All
/// MeatAxe applications should use this framework to achieve a consistent
/// behaviour. The functions described above enable an application to process
/// command line options and arguments. They also implement a rudimentary
/// on-line help and automatic processing of some common options like "-V".
/// Typically, a MeatAxe application's @c main() function should perform
/// the following steps:
/// - Create an application object using appAlloc().
/// - Process all command line options.
/// - Process simple options without arguments.
/// - Process options that take an argument. Options taking an integer
///   argument can be handeled with |appGetIntOption()|. Other option
///   arguments must be processed by the application.
/// - Call |appGetArgs()| to process the remaining command line arguments.
/// - Do whatever the application is supposed to do.
///
/// Here is a short example:
/// @code
/// MtxApplicationInfo_t AppInfo = { "sample", "MeatAxe sample application",
///     "\nSYNTAX\n"
///     "    sample [--all]] [--level <level>] <input> <ouput>\n"
/// };
/// MtxApplication_t *App;
/// int DoAll, Level;
/// char *Input, *Output;
///
/// int main(int argc, char *argv[])
/// {
///     App = appAlloc(&AppInfo,argc,argv);
///     DoAll = appGetOption(App,"-a --all");
///     Level = appGetIntOption(App,"-l --level",42,0,100);
///     appGetArguments(App,2,2);
///     Input = App->ArgV[0];
///     Output = App->ArgV[1];
///     DoIt();
///     appFree(App);
///     return 0;
/// }
/// @endcode
/// This sample application expects two arguments and recognizes two options.
/// One option ("-l") has an additional integer argument, @a level. The value
/// of @a level must be between 0 and 100. If not specified by the user,
/// a default value of 42 is used.

////////////////////////////////////////////////////////////////////////////////////////////////////

static int CheckForLongOption(MtxApplication_t *app, int i, const char *long_name)
{
   if (*long_name == 0) {
      return -1;
   }
   if (strcmp(app->origArgV[i] + 2,long_name)) {
      return -1;
   }
   app->isDone[i] = 0xFFFFFFFF;
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static int CheckForShortOption(MtxApplication_t *app, int i, char short_name, int needs_arg)
{
   const char *tab = app->origArgV[i] + 1;
   int k;
   for (k = 0; tab[k] != 0; ++k) {
      if (IS_DONE_1(app,i,k)) {
         continue;
      }
      if (tab[k] != short_name) {
         continue;
      }
      if (needs_arg && ((k > 0) || (tab[k + 1] != 0))) {
         mtxAbort(NULL,"Option '-%c' cannot be combined with other options", short_name);
         MARK_DONE(app,i);
         return -1;
      }
      MARK_DONE_1(app,i,k);
      return 0;
   }
   return -1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static int GetArg(MtxApplication_t *app, int i)
{
   if ((i >= app->optEnd - 1) || (app->isDone[i + 1] != 0)) {
      mtxAbort(NULL,"Missing argument after '%s'",app->origArgV[i]);
      return -1;
   }
   app->optArg = app->origArgV[i + 1];
   app->isDone[i + 1] = 0xFFFFFFFF;
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static int Find(MtxApplication_t *app, char short_name, const char *long_name,
                int needs_arg)
{
   int i;
   for (i = 0; i < app->optEnd; ++i) {
      int rc;

      if (IS_DONE(app,i)) {
         continue;
      }
      if (*app->origArgV[i] != '-') {
         continue;
      }
      if (app->origArgV[i][1] == '-') {
         rc = CheckForLongOption(app,i,long_name);
      } else {
         rc = CheckForShortOption(app,i,short_name,needs_arg);
      }
      if (rc == 0) {
         if (needs_arg) {
            if (GetArg(app,i) != 0) {
               return -1;
            }
         }
         app->optInd = i;
         return 0;
      }
   }
   return -1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Find command line option

static int FindSpec(MtxApplication_t *app, const char *spec, int needs_arg)
{
   const char *c;
   const char *short_name = "", *long_name = "";
   const char *err_text = "Invalid option specification";

   for (c = spec; *c != 0 && isspace((unsigned char)*c); ++c) {
   }
   if (*c != '-') {
      mtxAbort(MTX_HERE,"%s", err_text);
      return -1;
   }
   if (c[1] != '-') {
      short_name = c + 1;
      while (*c != 0 && !isspace((unsigned char)*c)) {
         ++c;
      }
      while (*c != 0 && isspace((unsigned char)*c)) {
         ++c;
      }
   }
   if (*c != 0) {
      if ((c[0] != '-') || (c[1] != '-')) {
         mtxAbort(MTX_HERE,"%s", err_text);
         return -1;
      }
      long_name = c + 2;
      while (*c != 0 && !isspace((unsigned char)*c)) {
         ++c;
      }
      if (*c != 0) {
         mtxAbort(MTX_HERE,"%s", err_text);
         return -1;
      }
   }

   return Find(app,*short_name,long_name,needs_arg);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Print the help text.

static void PrintHelp(const MtxApplicationInfo_t *ai)
{
   char v[100], *c;

   strcpy(v,mtxVersion());
   for (c = v; *c != 0 && *c != '$'; ++c) {
   }
   *c = 0;

   if (ai == NULL) {
      printf("MeatAxe Version %s\nNo help text available.\n",v);
   } else {
      printf("NAME\n    %s - %s\n    Version %s\n\n",
             ai->name,ai->description,v);
      printf("%s\n",ai->help);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Print the version info.

static void PrintVersion(const MtxApplicationInfo_t *ai)
{
   printf("%s\n",mtxVersion());
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Initialize the application.
/// This function initializes a MeatAxe application. It should be called
/// by the application's |main()| function before any other MeatAxe
/// library function is used. |argc| and |argv| should be the same values
/// that have been passed to |main()|. |ai|, if not NULL, must point to an
/// initialized application information structure.
///
/// %appAlloc() performs the following actions:
/// - It evaluates the MTXLIB environment variable. If set,
///   these variable overwrites the default library directory.
/// - It calls MtxInitLibrary() to initialize the MeatAxe library.
/// - It parses the command line and initializes internal data structures
///   for use with appGetOption() and related functions.
/// - Common options such as "-V" (see @ref prog_stdopts) are immediately evaluated.
///   These options are not visible to the application.
/// @param ai Application information.
/// @param argc Number of command line arguments.
/// @param argv List of command line arguments.
/// @return Pointer to the application object, or 0 on error.

MtxApplication_t *appAlloc(MtxApplicationInfo_t const *ai, int argc, char **argv)
{
   MtxApplication_t *a;
   const char *c;
   int time_limit;
   int i;

   a = ALLOC(MtxApplication_t);
   memset(a,0,sizeof(*a));
   a->context = mtxBegin(MTX_HERE, "Running program: %s", ai ? ai->name : "(no name)");

   /* Save the command line for later use.
      ------------------------------------ */
   a->optEnd = a->origArgC = argc - 1;
   a->origArgV = argv + 1;
   memset(a->isDone,0,sizeof(a->isDone));
   a->AppInfo = ai;

   /* Look for '--'.
      -------------- */
   for (i = 0; i < a->origArgC; ++i) {
      if (!strcmp(a->origArgV[i],"--")) {
         a->optEnd = i;
         a->isDone[i] = 0xFFFFFFFF;
         break;
      }
   }

   MtxMessageLevel = appGetCountedOption(a,"-V --verbose");
   MtxMessageLevel -= appGetCountedOption(a,"-Q --quiet");

   // Initialize the library
   if ((c = appGetTextOption(a,"-L --mtxlib",NULL)) != NULL && *c != 0) {
      setenv("MTXLIB", strdup(c), 1);
   }
   mtxInitLibrary(argv[0]);

   // Display help text
   if (appGetOption(a,"-h --help")) {
      PrintHelp(ai);
      exit(0);
   }
   if (appGetOption(a,"--version")) {
      PrintVersion(ai);
      exit(0);
   }
  
   // Handle common options.
   if ((time_limit = appGetIntOption(a,"-T --lime-limit",0,0,1000000)) > 0) {
      sysSetTimeLimit(time_limit);
   }
#if defined MTX_DEFAULT_THREADS
   const int nThreads = appGetIntOption(a,"-j --threads", MTX_DEFAULT_THREADS, 0, 1024);
   if (nThreads > 0)
      pexInit(nThreads);
#endif

   return a;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// End an application.
/// This function terminates a MeatAxe application. It should be called when the program ends.

void appFree(MtxApplication_t *a)
{
   long t = sysTimeUsed();
   MESSAGE(1,("%s: %ld.%ld seconds\n",a->AppInfo != NULL ?
              a->AppInfo->name : "meataxe",t / 10,t % 10));
   if (a->context > 0) mtxEnd(a->context);
   sysFree(a);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check for command line option.
/// This function checks if an option is present on the command line. Before
/// using appGetOption() you must call appAlloc() to initialize
/// the command line parser. Only simple options, i.e., options without an
/// argument can be recognized by this function. The return value is 1, if
/// the option is present and 0 otherwise.
///
/// The argument @a spec contains one or more names of the requested options.
/// If there is more than one name, names must be separated by spaces. All
/// names are considered equivalent, the user may use any of the names. The
/// leading "-" must always be included in the name. Typically an option has
/// only one or two names, as in the following example:
/// @code
/// int PrintAll = appGetOption(App,"-a --all");
/// @endcode
/// Here, the user may specify either '-a' or '--all' on the command line.
///
/// Each call of appGetOption() consumes at most one command line option.
/// If the user specifies the same option multiple times, only the first option
/// is recognized. Since this is usually not intended, an
/// appropriate error message will be generated when the application calls
/// appGetArguments(). If an option can be repeated on the command line
/// the application must make sure that all options are processed:
/// @code
/// while (appGetOption(App,"-a --all"))
/// { ... }
/// @endcode
/// The same remark applies to appGetIntOption() and appGetTextOption().
/// You may also use appGetCountedOption() to achieve a similar behaviour.
/// @param app Pointer to the application object.
/// @param spec The option name(s), see below.
/// @return 1 if the option is present, 0 otherwise.

int appGetOption(MtxApplication_t *app, const char *spec)
{
   return FindSpec(app,spec,0) == 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Count repeatable command line option.
/// This function counts how often a specific option is present on the
/// command line. For example, if the command line is
/// @code
/// sample -a -a -aa input output
/// @endcode
/// then MtxGetCountedOption("-a") returns 4. As with all command line
/// processing functions, you must call appAlloc() before
/// using this function.
///
/// Note: This function is included for compatibility reasons only.
/// New applications should not use it.
/// @param app Pointer to the application object.
/// @param spec The option name(s), see below.
/// @return Number of times the option is present on the command line.

int appGetCountedOption(MtxApplication_t *app, const char *spec)
{
   int count = 0;

   while (FindSpec(app,spec,0) == 0) {
      ++count;
   }
   return count;
}


/// Check for command line option.
/// This function checks if an option is present on the command line. The
/// option is expected to have a text value, and a pointer to the value is
/// returned. On the command line, the value must be separated from the
/// option by one or more spaces.
/// If the option is present on the command line but has no value, an
/// appropriate error message is generated. If the option is not present
/// on the command line, the function returns @a dflt as a default value.
/// @param app Pointer to the application object.
/// @param spec A list of names for the option, separated by spaces. See appGetOption().
/// @param dflt Default value.
/// @return Value of the option, or @a dflt if the option is not present.

const char *appGetTextOption(MtxApplication_t *app, const char *spec, const char *dflt)
{
   if (FindSpec(app,spec,1) != 0) {
      return dflt;
   }
   return app->optArg;
}


static int IsInteger(const char *c)
{
   if (*c == '-') {
      ++c;
   }
   if (!isdigit(*c)) {
      return 0;
   }
   do {
      ++c;
   } while (isdigit(*c));
   if (*c != 0) {
      return 0;
   }
   return 1;
}


/// Check for integer-valued option.
/// This function checks if an option is present on the command line. The
/// option is expected to have an integer value, and the value is returned.
/// On the command line, the value must be separated from the
/// option by one or more spaces.
/// If the option is present on the command line but has no value, an
/// appropriate error message is generated. If the option is not present
/// on the command line, the function returns @a dflt as a default value.
///
/// If the value on the command line is not within the range defined by @a min
/// and  @a max, an error message is generated. However, if @a min is
/// greater than @a max, no range check is performed.
/// @param app Pointer to the application object.
/// @param spec A list of names for the option, separated by spaces. See appGetOption().
/// @param dflt Default value.
/// @param min Minimum value of the option.
/// @param max Maximum value of the option.
/// @return Value of the option, or @a dflt if the option is not present.

int appGetIntOption(MtxApplication_t *app, const char *spec, int dflt,
                    int min, int max)
{
   const char *txt;
   int i;

   txt = appGetTextOption(app,spec,NULL);
   if (txt == NULL) {
      return dflt;
   }
   if (!IsInteger(txt)) {
      mtxAbort(NULL,"Invalid number after '%s'",app->origArgV[app->optInd]);
      return dflt;
   }
   i = atoi(txt);
   if ((min <= max) && ((i < min) || (i > max))) {
      mtxAbort(NULL,"Value after '%s' is out of range (%d..%d)",
                 app->origArgV[app->optInd],min,max);
      return dflt;
   }
   return i;
}


static int CheckDone(MtxApplication_t *app, int i)
{
   const char *c = app->origArgV[i];

   if (app->isDone[i] == 0xFFFFFFFF) {
      return 0;
   }
   if (*c != '-') {
      return 1;                 /* 1 = Ende der Optionen */
   }
   if (c[1] == '-') {
      if (app->isDone[i] != 0xFFFFFFFF) {
         mtxAbort(NULL,"Unknown option '%s', try --help",c);
         return -1;             /* -1 = Nicht ausgewertete Option */
      }
   } else {
      int k;
      ++c;
      for (k = 0; c[k] != 0; ++k) {
         if (!IS_DONE_1(app,i,k)) {
            char tmp[2];
            tmp[0] = c[k];
            tmp[1] = 0;
            mtxAbort(NULL,"Unknown option '-%s', try --help",tmp);
            return -1;                  /* -1 = Nicht ausgewertete Option */
         }
      }
   }
   return 0;                    /* 0 = Ok */
}


/// Get command line arguments.
/// This function must be called after all command line options have been
/// processed (see, for example, appGetOption()). The remaining words on
/// the command line are treated as non-optional arguments, the global variable
/// @c ArgV is set to the first argument and the number of arguments is
/// stored in @c |ArgC. An error message is generated, if there are unprocessed
/// options on the command line, or if the number of arguments is outside the
/// range specified by @a min_argc and @a max_argc.
/// @param app Pointer to the application object.
/// @param min_argc Minimum number of arguments expected.
/// @param max_argc Maximum number of arguments expected.
/// @return Number of arguments, or -1 on error.

int appGetArguments(MtxApplication_t *app, int min_argc, int max_argc)
{
   int i;

   // Check for unprocessed options
   for (i = 0; i < app->optEnd; ++i) {
      int rc = CheckDone(app,i);
      if (rc < 0) {
         return -1;
      }
      if (rc == 1) {
         break;
      }
   }

   // handle "--"
   if ((i == app->optEnd) && (app->optEnd < app->origArgC)) {
      ++i;
   }

   app->argC = app->origArgC - i;
   app->argV = app->origArgV + i;

   // check for options in argument list
   for (++i; i < app->origArgC; ++i) {
      if (app->isDone[i] != 0) {
         mtxAbort(MTX_HERE,"Option '%s' following non-optional argument", app->origArgV[i]);
         return -1;
      }
   }

   // check number of arguments
   if ((app->argC < min_argc) || (app->argC > max_argc)) {
      mtxAbort(MTX_HERE,"Invalid number of arguments, try --help");
      return -1;
   }
   return app->argC;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
