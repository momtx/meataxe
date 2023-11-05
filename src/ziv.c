////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Invert a matrix or permutation.
////////////////////////////////////////////////////////////////////////////////////////////////////



#include <stdlib.h>
#include "meataxe.h"



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


static const char *iname, *oname;

static MtxApplicationInfo_t AppInfo = { 
"ziv", "Invert Matrix or Permutation",
"SYNTAX"
"    ziv [-QV] <File> <Result>\n"
"\n"
"ARGUMENTS\n"
"    <File> .................. Input file name\n"
"    <Result> ................ Output file name\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"\n"
"FILES\n"
"    <Matrix> ................ I A square matrix or permutation\n"
"    <Result> ................ O Inverse matrix or permutation\n"
};







static MtxApplication_t *App = NULL;





static int Init(int argc, char **argv)

{
    /* Process command line options and arguments.
       ------------------------------------------- */
    if ((App = appAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (appGetArguments(App,2,2) != 2)
	return -1;
    iname = App->argV[0];
    oname = App->argV[1];
    return 0;
}



static void Cleanup()

{
    if (App != NULL)
	appFree(App);
}





int main(int argc, char **argv)
{
    void *x, *y;
    int rc = 0;

    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return -1;
    }
    x = objLoad(iname);
    if (x == NULL)
	return -1;
    y = objInverse(x);
    if (y != NULL)
    rc = objSave(y,oname);
    objFree(x);
    objFree(y);
    Cleanup();

    return rc;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_ziv ziv - Invert

@section ziv_syntax Command Line
<pre>
ziv @em Options @em Input @em Result
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par @em Input
  Input file.
@par @em Result
  Result file.

@section ziv_inp Input Files
@par @em Input
  A square natrix or permutation.

@section ziv_out Output Files
@par @em Result
  Inverted matrix or permutation.

@section ziv_desc Description
This program inverts a matrix or permutation.
**/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
