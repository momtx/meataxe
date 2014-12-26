////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Command line handling and application utilities
//
// (C) Copyright 1998-2014 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Global data

#if !defined(MTXBIN)
#define MTXBIN "/usr/local/mtx/bin"
#endif
#if !defined(MTXLIB)
#define MTXLIB "/usr/local/mtx/lib"
#endif

#define MTX_MAX_ARGS 150

   
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Local data

MTX_DEFINE_FILE_INFO
 
#define IS_DONE(app,i) (app->IsDone[i] == 0xFFFFFFFF)
#define IS_DONE_1(app,i,k) ((app->IsDone[i] & (1 << (k))) != 0)
#define MARK_DONE(app,i) (app->IsDone[i] = 0xFFFFFFFF)
#define MARK_DONE_1(app,i,k) (app->IsDone[i] |= (1 << (k)))


/// @defgroup app Application Framework
/// @{
/// The MeatAxe library provides a minimal framework for applications. All
/// MeatAxe applications should use this framework to achieve a consistent
/// behaviour. The functions described above enable an application to process 
/// command line options and arguments. They also implement a rudimentary
/// on-line help and automatic processing of some common options like "-V".
/// Typically, a MeatAxe application's @c main() function should perform 
/// the following steps:
/// - Create an application object using AppAlloc().
/// - Process all command line options.
/// - Process simple options without arguments.
/// - Process options that take an argument. Options taking an integer
///   argument can be handeled with |AppGetIntOption()|. Other option
///   arguments must be processed by the application.
/// - Call |AppGetArgs()| to process the remaining command line arguments.
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
///     App = AppAlloc(&AppInfo,argc,argv);
///     DoAll = AppGetOption(App,"-a --all");
///     Level = AppGetIntOption(App,"-l --level",42,0,100);
///     AppGetArguments(App,2,2);
///     Input = App->ArgV[0];
///     Output = App->ArgV[1];
///     DoIt();
///     AppFree(App);
///     return 0;
/// }
/// @endcode
/// This sample application expects two arguments and recognizes two options.
/// One option ("-l") has an additional integer argument, @a level. The value 
/// of @a level must be between 0 and 100. If not specified by the user,
/// a default value of 42 is used.
///
/// @par Temporary Directories  TODO: move to AppCreateTempDir()
/// If an application needs a temporary directory to store interdediate files,
/// it can use AppCreateTempDir(). This function creates a new directory and
/// returns its name. The directory will be removed when the application object
/// is destroyed, i.e., in |AppFree()|. Note that the application must not leave
/// any file in the temporary directory. Otherwise, a run-time error is generated
/// when AppFree() tries to remove the directory.
/// The following code example shows how to use AppCreateTempDir():
/// @code
/// int main(int argc, char **argv)
/// {
///     MtxApplication_t *app = AppAlloc(argc,argv,&appinfo);
///     const char *tmpdir = AppCreateTempDir(app);
///     ...
///
///     sprintf(file_name,"%s/%s",tempdir,"test");
///     MatSave(mat,fn);
///     ...
///     SysRemoveFile(fn);
///
///     AppFree(app);
///     return 0;
/// }
/// @endcode

////////////////////////////////////////////////////////////////////////////////////////////////////
/// MeatAxe Binary Directory.
/// This variable contains the name of the MeatAxe binary directory. This directory
/// is used, for example, to find the @ref prog_maketab "maketab" program for automatic
/// arithmetic table generation. The value of MtxBinDir can be set on the command line
/// with the "-B" option. Otherwise, the value of the environment variable
/// MTXBIN is used. If neither "-B" nor MTXBIN are defined, the default
/// directory, which was selected when building the MeatAxe, is used.

char MtxBinDir[250] = MTXBIN;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// MeatAxe Library Directory.
/// This variable contains the name of the MeatAxe library directory.
/// Arithmetic table files are searched in this directory. The value of 
/// MtxLibDir can be set on the command line with the "-L" option. Otherwise, 
/// the value of the environment variable MTXLIB is used. If neither "-L" nor 
/// MTXLIB are defined, the default directory, which was selected when 
/// building the MeatAxe, is used.
/// @see  MtxBinDir 

char MtxLibDir[250] = MTXLIB;

////////////////////////////////////////////////////////////////////////////////////////////////////

static int CheckForLongOption(MtxApplication_t *app, int i, const char *long_name)
{
    if (*long_name == 0)
	return -1;
    if (strcmp(app->OrigArgV[i]+2,long_name))
        return -1;
    app->IsDone[i] = 0xFFFFFFFF;
    return 0;
}
 
////////////////////////////////////////////////////////////////////////////////////////////////////
 
static int CheckForShortOption(MtxApplication_t *app, int i, char short_name, int needs_arg)
{
    const char *tab = app->OrigArgV[i] + 1;
    int k;
    for (k = 0; tab[k] != 0; ++k)
    {
        if (IS_DONE_1(app,i,k))
            continue;
        if (tab[k] != short_name)
            continue;
        if (needs_arg && (k > 0 || tab[k+1] != 0))
        {
	    MTX_ERROR1("Option '-%c' cannot be combined with other options",
	    	short_name);
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
    if (i >= app->OptEnd - 1 || app->IsDone[i+1] != 0)
    {
        MTX_ERROR1("Missing argument after '%s'",app->OrigArgV[i]);
        return -1;
    }
    app->OptArg = app->OrigArgV[i+1];
    app->IsDone[i+1] = 0xFFFFFFFF;
    return 0;
}
 
////////////////////////////////////////////////////////////////////////////////////////////////////

static int Find(MtxApplication_t *app, char short_name, const char *long_name, 
    int needs_arg)
 
{
    int i;
    for (i = 0; i < app->OptEnd; ++i)
    {
        int rc;
 
        if (IS_DONE(app,i))
            continue;
        if (*app->OrigArgV[i] != '-')
            continue;
        if (app->OrigArgV[i][1] == '-')
            rc = CheckForLongOption(app,i,long_name);
        else
            rc = CheckForShortOption(app,i,short_name,needs_arg);
        if (rc == 0)
        {
            if (needs_arg)
            {
                if (GetArg(app,i) != 0)
                    return -1;
            }
            app->OptInd = i;
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
 

    for (c = spec; *c != 0 && isspace((unsigned char)*c); ++c)
    	;
    if (*c != '-')
    {
        MTX_ERROR(err_text);
        return -1;
    }
    if (c[1] != '-')
    {
    	short_name = c + 1;
	while (*c != 0 && !isspace((unsigned char)*c))
	    ++c;
	while (*c != 0 && isspace((unsigned char)*c))
	    ++c;
    }
    if (*c != 0)
    {
    	if (c[0] != '-' || c[1] != '-')
    	{
            MTX_ERROR(err_text);
            return -1;
        }
	long_name = c + 2;
	while (*c != 0 && !isspace((unsigned char)*c))
	    ++c;
	if (*c != 0)
    	{
            MTX_ERROR(err_text);
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

    strcpy(v,MtxVersion);
    for (c = v; *c != 0 && *c != '$'; ++c);
    *c = 0;

    if (ai == NULL)
    {
	printf("MeatAxe Version %s\nNo help text available.\n",v);
    }
    else
    {
	printf("NAME\n    %s - %s\n    Revision %s\n\n",
	      ai->Name,ai->Description,v);
	printf("%s\n",ai->Help);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Initialize the application.
/// This function initializes a MeatAxe application. It should be called
/// by the application's |main()| function before any other MeatAxe
/// library function is used. |argc| and |argv| should be the same values
/// that have been passed to |main()|. |ai|, if not NULL, must point to an
/// initialized application information structure.
/// 
/// %AppAlloc() performs the following actions:
/// - It evaluates the MTXLIB and MTXBIN environment variables. If set,
///   these variables overwrite the default directories.
/// - It calls MtxInitLibrary() to initialize the MeatAxe library.
/// - It parses the command line and initializes internal data structures
///   for use with AppGetOption() and related functions.
/// - Common options such as "-V" (see @ref prog_stdopts) are immediately evaluated.
///   These options are not visible to the application.
/// @param ai Application information.
/// @param argc Number of command line arguments.
/// @param argv List of command line arguments.
/// @return Pointer to the application object, or 0 on error.

MtxApplication_t *AppAlloc(MtxApplicationInfo_t const *ai, int argc, const char **argv)
{
    MtxApplication_t *a;
    const char *c;
    int time_limit;
    int i;


    if ((a = ALLOC(MtxApplication_t)) == NULL)
	return NULL;
    memset(a,0,sizeof(*a));

    /* Save the command line for later use.
       ------------------------------------ */
    a->OptEnd = a->OrigArgC = argc - 1;
    a->OrigArgV = argv + 1;
    memset(a->IsDone,0,sizeof(a->IsDone));
    a->AppInfo = ai;

    /* Look for '--'.
       -------------- */
    for (i = 0; i < a->OrigArgC; ++i)
    {
	if (!strcmp(a->OrigArgV[i],"--"))
	{
	    a->OptEnd = i;
	    a->IsDone[i] = 0xFFFFFFFF;
	    break;
	}
    }

    /* Process environment variables
       ----------------------------- */
    if ((c = getenv("MTXBIN")) != NULL)
	strcpy(MtxBinDir,c);
    if ((c = getenv("MTXLIB")) != NULL)
	strcpy(MtxLibDir,c);

    /* Initialize the library
       ---------------------- */
    MtxInitLibrary();

    /* Display help text
       ----------------- */
    if (AppGetOption(a,"-h --help"))
    {
	PrintHelp(ai);
	exit(0);
    }
    
    /* Check for common options
       ------------------------ */
    MtxMessageLevel = AppGetCountedOption(a,"-V --verbose");
    if (AppGetOption(a,"-Q --quiet"))
	MtxMessageLevel = -1000;
    if ((c = AppGetTextOption(a,"-L --mtxlib",NULL)) != NULL)
	strcpy(MtxLibDir,c);
    if ((c = AppGetTextOption(a,"-B --mtxbin",NULL)) != NULL)
	strcpy(MtxBinDir,c);
    MtxOpt_UseOldWordGenerator = AppGetOption(a,"--old-word-generator");
    if ((time_limit = AppGetIntOption(a,"-T --lime-limit",0,0,1000000)) > 0)
	SysSetTimeLimit(time_limit);

    return a;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// End an application.
/// This function terminates a MeatAxe application. It should be called
/// immediately before the |main()| function exits. |AppFree()| removes
/// any temporary directory created with |AppCreateTempDir()|. If the
/// directory is not empty at this point, a run-time error will result.
/// @return 0 on success, -1 on error.

int AppFree(MtxApplication_t *a)
{
    long t = SysTimeUsed();

    /* Remove the temporary directory.
       ------------------------------- */
    if (a->TempDirName[0] != 0)
	SysRemoveDirectory(a->TempDirName);

    MESSAGE(1,("%s: %ld.%ld seconds\n",a->AppInfo != NULL ?
	a->AppInfo->Name : "meataxe",t/10,t%10));
    MtxCleanupLibrary();
    SysFree(a);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check for command line option.
/// This function checks if an option is present on the command line. Before
/// using AppGetOption() you must call AppAlloc() to initialize
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
/// int PrintAll = AppGetOption(App,"-a --all");
/// @endcode
/// Here, the user may specify either '-a' or '--all' on the command line.
///
/// Each call of AppGetOption() consumes at most one command line option.
/// If the user specifies the same option multiple times, only the first option
/// is recognized. Since this is usually not intended, an
/// appropriate error message will be generated when the application calls 
/// AppGetArguments(). If an option can be repeated on the command line
/// the application must make sure that all options are processed:
/// @code
/// while (AppGetOption(App,"-a --all"))
/// { ... }
/// @endcode
/// The same remark applies to AppGetIntOption() and AppGetTextOption().
/// You may also use AppGetCountedOption() to achieve a similar behaviour.
/// @param app Pointer to the application object.
/// @param spec The option name(s), see below.
/// @return 1 if the option is present, 0 otherwise.

int AppGetOption(MtxApplication_t *app, const char *spec)
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
/// processing functions, you must call AppAlloc() before
/// using this function.
///
/// Note: This function is included for compatibility reasons only.
/// New applications should not use it.
/// @param app Pointer to the application object.
/// @param spec The option name(s), see below.
/// @return Number of times the option is present on the command line.

int AppGetCountedOption(MtxApplication_t *app, const char *spec)
{
    int count = 0;

    while (FindSpec(app,spec,0) == 0)
	++count;
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
/// @param spec A list of names for the option, separated by spaces. See AppGetOption().
/// @param dflt Default value.
/// @return Value of the option, or @a dflt if the option is not present.

const char *AppGetTextOption(MtxApplication_t *app, const char *spec, const char *dflt)
{
    if (FindSpec(app,spec,1) != 0)
        return dflt;
    return app->OptArg;
}



static int IsInteger(const char *c)
{
    if (*c == '-') 
	++c;
    if (!isdigit(*c))
	return 0;
    do ++c; while (isdigit(*c));
    if (*c != 0)
	return 0;
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
/// @param spec A list of names for the option, separated by spaces. See AppGetOption().
/// @param dflt Default value.
/// @param min Minimum value of the option.
/// @param max Maximum value of the option.
/// @return Value of the option, or @a dflt if the option is not present.


int AppGetIntOption(MtxApplication_t *app, const char *spec, int dflt, 
    int min, int max)

{
    const char *txt;
    int i;

    txt = AppGetTextOption(app,spec,NULL);
    if (txt == NULL)
	return dflt;
    if (!IsInteger(txt))
    {
	MTX_ERROR1("Invalid number after '%s'",app->OrigArgV[app->OptInd]);
	return dflt;
    }
    i = atoi(txt);
    if (min <= max && (i < min || i > max))
    {
	MTX_ERROR3("Value after '%s' is out of range (%d..%d)",
	    app->OrigArgV[app->OptInd],min,max);
	return dflt;
    }
    return i;
}



static int CheckDone(MtxApplication_t *app, int i)

{
    const char *c = app->OrigArgV[i];

    if (app->IsDone[i] == 0xFFFFFFFF)
    	return 0;
    if (*c != '-')
    	return 1;		/* 1 = Ende der Optionen */
    if (c[1] == '-')
    {
    	if (app->IsDone[i] != 0xFFFFFFFF)
	{
	    MTX_ERROR1("Unknown option '%s', try --help",c);
	    return -1;		/* -1 = Nicht ausgewertete Option */
	}
    }
    else
    {
    	int k;
	++c;
	for (k = 0; c[k] != 0; ++k)
	{
	    if (!IS_DONE_1(app,i,k))
	    {
		char tmp[2];
		tmp[0] = c[k];
		tmp[1] = 0;
		MTX_ERROR1("Unknown option '-%s', try --help",tmp);
		return -1;		/* -1 = Nicht ausgewertete Option */
	    }
	}
    }
    return 0;			/* 0 = Ok */
}




/// Get command line arguments.
/// This function must be called after all command line options have been 
/// processed (see, for example, AppGetOption()). The remaining words on 
/// the command line are treated as non-optional arguments, the global variable 
/// @c ArgV is set to the first argument and the number of arguments is 
/// stored in @c |ArgC. An error message is generated, if there are unprocessed 
/// options on the command line, or if the number of arguments is outside the 
/// range specified by @a min_argc and @a max_argc.
/// @param app Pointer to the application object.
/// @param min_argc Minimum number of arguments expected.
/// @param max_argc Maximum number of arguments expected.
/// @return Number of arguments, or -1 on error.

int AppGetArguments(MtxApplication_t *app, int min_argc, int max_argc)
{
    int i;

    /* Pruefen, ob alle Optionen abgearbeitet wurden.
       ---------------------------------------------- */
    for (i = 0; i < app->OptEnd; ++i)
    {
    	int rc = CheckDone(app,i);
	if (rc < 0)
	    return -1;
	if (rc == 1)
	    break;
    }

    /* '--' ueberspringen, falls vorhanden.
       ------------------------------------ */
    if (i == app->OptEnd && app->OptEnd < app->OrigArgC)
    	++i;
    
    app->ArgC = app->OrigArgC - i;
    app->ArgV = app->OrigArgV + i;

    /* Pruefe, ob weitere Optionen folgen.
       ----------------------------------- */
    for (++i; i < app->OrigArgC; ++i)
    {
	if (app->IsDone[i] != 0)
	{
	    MTX_ERROR1("Option '%s' following non-optional argument",
		app->OrigArgV[i]);
	    return -1;
	}
    }

    /* Pruefe die Anzahl der Argumente.
       -------------------------------- */
    if (app->ArgC < min_argc || app->ArgC > max_argc)
    {
	MTX_ERROR("Invalid number of arguments, try --help");
	return -1;
    }
    return app->ArgC;
}


/// Create a temporary directory.
/// This function creates a new, empty directory for use by the application.
/// The location of this directory is system dependent. When the application
/// objects is destroyed, the temporary directory will be removed 
/// automatically, if there are no files left in it.
/// After successfull creation of the directory, the function returns the
/// directory name. The return value is a pointer to an internal buffer,
/// which must not be modified by the application. If the directory cannot be
/// created for some reason, the return value is 0.
///
/// If the application calls %AppCreateTempDir() more than once, only the
/// first call actually creates a directory. The later calls just return the
/// name of the directory already created.
/// @param app Pointer to the application object.
/// @return Directory name or 0 on error.

const char *AppCreateTempDir(MtxApplication_t *app)
{
    /* Check if we have already created a temporary directory.
       ------------------------------------------------------- */
    if (app->TempDirName[0] != 0)
	return app->TempDirName;

    /* Make the directory.
       ------------------- */
    sprintf(app->TempDirName,"mtxtmp.%5.5d",SysGetPid());
    if (SysCreateDirectory(app->TempDirName) != 0)
    {
	MTX_ERROR("Cannot create temporary directory");
	app->TempDirName[0] = 0;
	return NULL;
    }

    return app->TempDirName;
}


/// @}
