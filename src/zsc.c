////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Spinup with script
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <string.h>
#include <stdlib.h>



/* --------------------------------------------------------------------------
   Global data
   -------------------------------------------------------------------------- */


static MtxApplicationInfo_t AppInfo = { 
"zsc", "Spin up with script",
"\n"
"SYNTAX\n"
"    zsc [<Options>] <Gen> <Seed> <Op> [<Out>]\n"
"\n"
"ARGUMENTS\n"
"    <Gen> ................... Generator name\n"
"    <Seed> .................. Seed vector(s)\n"
"    <Op> .................... Spin-up script\n"
"    <Out> ................... Result\n"        
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -g <#Gen> ............... Set number of generators.\n"
"\n"
"FILES\n"
"    <Gen>.{1,2...} .......... I  Generators\n"
"    <Seed> .................. I  Seed vectors\n"
"    <Op> .................... I  Spin-up script\n"
"    <Out> ................... I  Output file\n"
};

static MtxApplication_t *App = NULL;

static const char *GenName = NULL;	/* File name of generators */
static const char *SeedName = NULL;	/* Seed vector file name */
static const char *OutName = NULL;	/* Output file name (base name) */
static const char *OpName = NULL;	/* Name of <op> file */

static int fl, nor;		/* Field and size of the generators */
static IntMatrix_t *OpTable;	/* The spin-up script */
static int ngen = 2;		/* Number of generators */
static MatRep_t *Rep;		/* The generators */
static Matrix_t *Seed;		/* Seed vectors */

#define OPVEC(i) OpTable->data[i * 2]
#define OPGEN(i) OpTable->data[i * 2 + 1]

/* --------------------------------------------------------------------------
   ReadFiles() - Read input files

   Description:
     This function reads the generators and the script file. It also
     allocates a workspace for spin-up.
   -------------------------------------------------------------------------- */
	     
static int ReadFiles()

{
    int i;

    /* Read the generators.
      --------------------- */
    Rep = mrLoad(GenName,ngen);
    if (Rep == NULL)
	return -1;
    fl = Rep->Gen[0]->field;
    nor = Rep->Gen[0]->nor;

    /* Read the script.
       ---------------- */
    MTX_LOGD("Reading %s",OpName);
    OpTable = imatLoad(OpName);
    if (OpTable == NULL)
	return -1;
    ConvertSpinUpScript(OpTable);   /* 2.3 compatibility */
    if (OpTable->noc != 2)
    {
	mtxAbort(MTX_HERE,"%s: bad number of columns",OpName);
	return -1;
    }
   /* Das ist doch unsinnig! J.M. */
   /* if (OpTable->nor != nor)
    {
	mtxAbort(MTX_HERE,"%s: bad number of rows",OpName);
	return -1;
    } */

    /* Check the script for errors
       --------------------------- */
    /* if (OPVEC(0) != 0 || OPGEN(0) != -1) nanana, J.M. */
    if (OPGEN(0) != -1)
	mtxAbort(MTX_HERE,"Illegal script (does not start with seed vector)");
    if (OPVEC(0) != 1)
	MTX_LOGI("Note: script does not start with first vector");
    for (i = 1; i < OpTable->nor; ++i)
    {
	if (OPGEN(i) == -1)
	    mtxAbort(MTX_HERE,"Illegal skript (more than 1 seed vector)");
	if (OPGEN(i) < 0 || OPGEN(i) >= ngen)
	    mtxAbort(MTX_HERE,"Illegal skript (pos %d: generator out of range)",i);
	/* if (OPVEC(i) <= 0 || OPVEC(i) >= i) nanana, J.M. */
	if (OPVEC(i) < 0 || OPVEC(i) >= i)
	    mtxAbort(MTX_HERE,"Illegal skript (pos %d: vector out of range)",i);
    }

    /* Read the seed vectors.
       ---------------------- */
    Seed = matLoad(SeedName);
    if (Seed == NULL)
	return -1;
    if (Seed->noc != nor || Seed->field != fl)
    {
	mtxAbort(MTX_HERE,"%s.1 and %s: %s",GenName,SeedName,MTX_ERR_INCOMPAT);
	return 1;
    }
    MTX_LOGD("%s: %d seed vectors",SeedName,Seed->nor);
    return 0;
}



/* ------------------------------------------------------------------
   Init() - Program initialization
   ------------------------------------------------------------------ */

static int Init(int argc, char **argv)

{
    if ((App = appAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;

    /* Options.
       -------- */
    ngen = appGetIntOption(App,"-g",2,1,1000);

    /* Arguments.
       ---------- */
    if (appGetArguments(App,3,4) < 0)
	return -1;
    OpName = App->argV[2];
    SeedName = App->argV[1];
    GenName = App->argV[0];
    OutName = App->argC > 3 ? App->argV[3] : SeedName;

    /* Read input files.
       ----------------- */
    if (ReadFiles() != 0)
	return -1;

    return 0;
}






/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char *argv[])
{
    int seedno;

    if (Init(argc, argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return 1;
    }

    /* Spinup each seed vector.
       ------------------------ */
    for (seedno = 0; seedno < Seed->nor; ++seedno)
    {
	char fn[200];
	Matrix_t *seed_vec = matCutRows(Seed,seedno,1);
	Matrix_t *result = SpinUpWithScript(seed_vec,Rep,OpTable);
        if (result == NULL)
	    return 1;
	sprintf(fn,"%s.%d",OutName,seedno+1);
	matSave(result,fn);
	matFree(result);
	matFree(seed_vec);
    }

    /* Clean up.
       --------- */
    appFree(App);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_zsc zsc - Spin-Up With Script

@section zsc_syntax Command Line
<pre>
zsc [@em Options] [-g @em NGen] @em Gen @em Seed @em Script @em Output
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts.
@par -g @em NGen
  Set the number of generators (default: 2).
@par @em Gen
  Generator base name.
@par @em Seed
  Seed space.
@par @em Script
  Spin-up script.
@par @em Output
  Output base name. Default is @em Script.

@section zsc_inp Input Files
@par @em Gen.1, @em Gen.2, ...
  Generators.
@par @em Seed
  Seed space.
@par @em Op
  Spin-up script.

@section zsc_out Output Files
@par @em Output.0001, @em Output.0002, ...
  Spin-Up result, one file per seed vector.

@section zsc_desc Description
This program reads in two or more matrices (generators), a list of seed vectors 
and a list of operations (the script). Then, ZSC applies the script to each 
seed vector and writes the output in a separate file for each vector.

The generators must be square matrices over the same field. @em Seed 
must be a matrix over the same field with the same number of columns. 
@em Script must be an integer matrix with two columns in the format produced
by the @ref prog_zsp "zsp" program. 
Only one seed vectur may be used in the script, i.e., each row except the
first row of @em Script must be of the form (x,y) with x≥0.

The number of generators is 2 by default. This can be changed
by using the "-g" option. The number of output files equals
the number of seed vectors. If no fourth argument is given,
the output name defaults to @em Seed. For example,
<pre>zsc -g 3 gen seed op</pre>
reads three generators fom "gen.1", "gen.2" and "gen.3",
seed vectors from `seed', the script from "op", and
writes the output to "seed.0001", "seed.0002", ...
*/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
