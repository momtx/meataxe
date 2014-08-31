/* ============================= C MeatAxe ==================================
   File:        $Id: zpo.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Power of a matrix or permutation.
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"

#include <ctype.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

static int Power;

static MtxApplicationInfo_t AppInfo = { 
"power", "Power of a matrix or permutation", 
"SYNTAX\n"
"    zpo [<Options>] <In> <n> <Result>\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"\n"
"ARGUMENTS\n"
"    <In> ........... Input: Matrix or permutation\n"
"    <n> ............ Power to compute (e.g., 5 or pwr5)\n"
"    <Result> ....... Ouput: <n>-th power of <In>\n"
};

static MtxApplication_t *App = NULL;








/* ------------------------------------------------------------------
   init() - Process command line options and arguments
   ------------------------------------------------------------------ */

static int Init(int argc, const char **argv)

{
    /* Process command line options.
       ----------------------------- */
    App = AppAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;

    /* Process arguments.
       ------------------ */
    if (AppGetArguments(App,3,3) < 0)
	return -1;
    if (!strncmp(App->ArgV[1],"pwr",3))
        Power = atoi(App->ArgV[1] + 3);		/* ZSM compatibility (pwrN) */
    else
        Power = atoi(App->ArgV[1]);

    return 0;
}




static void Cleanup()

{
    AppFree(App);
}



static int CalcPower()
{
    void *x = XLoad(App->ArgV[0]);
    void *y;

    if (x == NULL)
	return -1;
    y = XPower(x,Power);
    if (y == NULL)
	return -1;
    XSave(y,App->ArgV[2]);
    XFree(x);
    XFree(y);
    return 0;

}






/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, const char **argv)
{
    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return 1;
    }
    CalcPower();
    Cleanup();
    return EXIT_OK;
}


/**
@page prog_zpo zpo - Power

<<<<<<< HEAD
@section syntax Command Line
=======
@section zpo_syntax Command Line
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
<pre>
zpo [@em Options] @em Input @em N @em Result
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par @em Input
  Input file.
@par @em N
  Exponent.
@par @em Result
  Output file.

<<<<<<< HEAD
@section inp Input Files
@par @em Input
  Input file.

@section inp Output Files
@par @em Result
  Output file.

@section desc Description
=======
@section zpo_inp Input Files
@par @em Input
  Input file.

@section zpo_out Output Files
@par @em Result
  Output file.

@section zpo_desc Description
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
This program calculates the @em N-th power of a matrix or permutation.
The exponent @em N may be negative. For compatibility with the old @b zsm program,
the power may be specified in the @b zsm format. The following example
shows two equivalent ways to calculate the 69-th power of a matrix:
<pre>
zpo matrix 69 result
zpo matrix pwr69 result
</pre>
**/

