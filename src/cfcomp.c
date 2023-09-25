////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - cfcomp: Compare irreducible constituents
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>
#include <stdlib.h>

static MtxApplicationInfo_t AppInfo = {
   "cfcomp", "Compare irreducible constituents",
   "SYNTAX\n"
   "    cfcomp <Module> <Module2> ...\n"
   "\n"
   "ARGUMENTS\n"
   "    <Module> ................ The reference module (must be chopped).\n"
   "    <ModuleN> ............... Another irreducible module (only generators needed).\n"
   "\n"
   "OPTIONS\n"
   MTX_COMMON_OPTIONS_DESCRIPTION
   "\n"
   "FILES\n"
   "    <Module>.cfinfo.......... I  Constituent information (generated by CHOP)\n"
};

static MtxApplication_t *app = NULL;
static Lat_Info infoA;
static MatRep_t *irredA[LAT_MAXCF];
static MatRep_t *repB;

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
   char fn[100];
   int i;

   app = appAlloc(&AppInfo,argc,argv);
   appGetArguments(app,2,2000);
   latReadInfo(&infoA,app->argV[0]);

   // Read the generators for each composition factor
   for (i = 0; i < infoA.nCf; ++i) {
      sprintf(fn,"%s%s",infoA.BaseName,latCfName(&infoA,i));
      MESSAGE(1,("Reading %s\n",fn));
      irredA[i] = mrLoad(fn,infoA.NGen);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void cleanup()
{
   for (int i = 0; i < infoA.nCf; ++i) {
      mrFree(irredA[i]);
   }
   appFree(app);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void findEquiv(const char *name)
{
   for (int i = 0; i < infoA.nCf; ++i) {
      if (repB->Gen[0]->nor != irredA[i]->Gen[0]->nor) {
         continue;
      }
      if (IsIsomorphic(irredA[i],infoA.Cf + i,repB,NULL,0)) {
         MESSAGE(0,("%s = %s%s\n",name,infoA.BaseName,latCfName(&infoA,i)));
         return;
      }
   }
   MESSAGE(0,("%s not found in %s\n",name,infoA.BaseName));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void compare(const char *nameB)
{
   repB = mrLoad(nameB,infoA.NGen);
   findEquiv(nameB);
   mrFree(repB);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
   init(argc,argv);
   for (int i = 1; i < app->argC; ++i) {
      compare(app->argV[i]);
   }
   cleanup();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// @page prog_cfcomp cfcomp - Compare Irreducible Constituents
///
/// @section cfcomp_syntax Command Line
/// <pre>
/// cfcomp [@em Options] @em Module @em Irred [@em Irred ...]
/// </pre>
///
/// @par @em Options
/// Standard options, see @ref prog_stdopts
/// @par @em Module
/// Name of the representation.
/// @par @em Irred
/// Irreducible module.
///
/// @section cfcomp_inp Input Files
/// @par @em Module.cfinfo
/// Constituent info file.
/// @par @em Irred.1, @em Irred.2, ...
/// Generators.
///
/// @section cfcomp_desc Description
/// After @em Module has been chopped, you can use this program to determine
/// if a given irreducible module, @em Irred, is a constituent of @em Module.
/// If yes, the program finds out which of the constituents of @em Module is
/// isomorphic to @em Irred.
///
/// The program needs at least two arguments. The first argument is the name of
/// the chopped module. The remaining arguments are treated as names of irreducible
/// modules. Each of these irreducible modules is checked against the irreducible
/// constituents of @em Module.

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
