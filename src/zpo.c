////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Power of a matrix or permutation.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <ctype.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


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

static int Init(int argc, char **argv)

{
    /* Process command line options.
       ----------------------------- */
    App = appAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;

    /* Process arguments.
       ------------------ */
    if (appGetArguments(App,3,3) < 0)
	return -1;
    if (!strncmp(App->argV[1],"pwr",3))
        Power = atoi(App->argV[1] + 3);		/* ZSM compatibility (pwrN) */
    else
        Power = atoi(App->argV[1]);

    return 0;
}




static void Cleanup()

{
    appFree(App);
}



static int CalcPower()
{
    void *x = objLoad(App->argV[0]);
    void *y;

    if (x == NULL)
	return -1;
    y = objPower(x,Power);
    if (y == NULL)
	return -1;
    objSave(y,App->argV[2]);
    objFree(x);
    objFree(y);
    return 0;

}






/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return 1;
    }
    CalcPower();
    Cleanup();
    return EXIT_OK;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_zpo zpo - Power

@section zpo_syntax Command Line
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

@section zpo_inp Input Files
@par @em Input
  Input file.

@section zpo_out Output Files
@par @em Result
  Output file.

@section zpo_desc Description
This program calculates the @em N-th power of a matrix or permutation.
The exponent @em N may be negative. For compatibility with the old @b zsm program,
the power may be specified in the @b zsm format. The following example
shows two equivalent ways to calculate the 69-th power of a matrix:
<pre>
zpo matrix 69 result
zpo matrix pwr69 result
</pre>
**/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
