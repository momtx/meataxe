////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Compare irreducibe constituents
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <meataxe.h>
#include <string.h>
#include <stdlib.h>


/* ------------------------------------------------------------------
   Global variables
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

static MtxApplicationInfo_t AppInfo = { 
"cfcomp", "Compare irreducible constituents",
"SYNTAX\n"
"    cfcomp <Module <Module2> ...\n"
"\n"
"ARGUMENTS\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"\n"
"FILES\n"
"    <Name>.cfinfo............ I  Constituent information (generated by CHOP)\n"
};

static MtxApplication_t *App = NULL;
static Lat_Info InfoA;
static MatRep_t *CfGenA[LAT_MAXCF];
static MatRep_t *GenB;


/* ------------------------------------------------------------------
   init ()
   ------------------------------------------------------------------ */

static int Init(int argc, const char **argv)

{
    char fn[100];
    int i;

    App = AppAlloc(&AppInfo,argc,argv);
    if (AppGetArguments(App,2,2000) < 0)
	return -1;
    if (Lat_ReadInfo(&InfoA,App->ArgV[0]) != 0)
	return -1;

    /* Read the generators for each composition factor
       ----------------------------------------------- */
    for (i = 0; i < InfoA.NCf; ++i)
    {
	sprintf(fn,"%s%s",InfoA.BaseName,Lat_CfName(&InfoA,i));
	MESSAGE(1,("Reading %s\n",fn));
    	if ((CfGenA[i] = MrLoad(fn,InfoA.NGen)) == NULL)
	    return -1;
    }
    return 0;
}



static void Cleanup()
{
    int i;
    for (i = 0; i < InfoA.NCf; ++i)
	MrFree(CfGenA[i]);
    AppFree(App);
}



void ReadGens(const char *name)

{
    GenB = MrLoad(name,InfoA.NGen);
}


void FreeGens()

{
    MrFree(GenB);
}






/* ------------------------------------------------------------------
   FindEquiv()
   ------------------------------------------------------------------ */

static void FindEquiv(const char *name)

{
    int i;

    for (i = 0; i < InfoA.NCf; ++i)
    {
	/* The dimension must be equal. */
	if (GenB->Gen[0]->Nor != CfGenA[i]->Gen[0]->Nor)
	    continue;

	/* Check for equivalence. */
	if (IsIsomorphic(CfGenA[i],InfoA.Cf + i,GenB,NULL,0))
	{
	    MESSAGE(0,("%s = %s%s\n",name,InfoA.BaseName,Lat_CfName(&InfoA,i)));
	    return;
	}
    }
    MESSAGE(0,("%s not found in %s\n",name,InfoA.BaseName));
}


/* ------------------------------------------------------------------
   Compare()
   ------------------------------------------------------------------ */

static void Compare(const char *name)
{
    ReadGens(name);
    FindEquiv(name);
    FreeGens(name);
}



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, const char *argv[])
{   
    int i;

    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return -1;
    }
    for (i = 1; i < App->ArgC; ++i)
	Compare(App->ArgV[i]);
    Cleanup();
    return 0;
}



/**
@page prog_cfcomp cfcomp - Compare Irreducible Constituents

@section cfcomp_syntax Command Line
<pre>
cfcomp [@em Options] @em Module @em Irred [@em Irred ...]
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par @em Module
  Name of the representation.
@par @em Irred
  Irreducible module.

@section cfcomp_inp Input Files
@par @em Module.cfinfo
  Constituent info file.
@par @em Irred.1, @em Irred.2, ...
  Generators.

@section cfcomp_desc Description
After @em Module has been chopped, you can use this program to determine 
if a given irreducible module, @em Irred, is a constituent of @em Module.
If yes, the program finds out which of the constituents of @em Module is 
isomorphic to @em Irred.

The program needs at least two arguments. The first argument is the name of
the chopped module. The remaining arguments are treated as names of irreducible
modules. Each of these irreducible modules is checked against the irreducible
constituents of @em Module.

**/
