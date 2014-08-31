/* ============================= C MeatAxe ==================================
   File:        $Id: ziv.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Invert a matrix or permutation.
   --------------------------------------------------------------------------
   (C) Copyright 1997 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */



#include <stdlib.h>
#include "meataxe.h"



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

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





static int Init(int argc, const char **argv)

{
    /* Process command line options and arguments.
       ------------------------------------------- */
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (AppGetArguments(App,2,2) != 2)
	return -1;
    iname = App->ArgV[0];
    oname = App->ArgV[1];
    return 0;
}



static void Cleanup()

{
    if (App != NULL)
	AppFree(App);
}





int main(int argc, const char **argv)

{
    void *x, *y;
    int rc = 0;

    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return -1;
    }
    x = XLoad(iname);
    if (x == NULL)
	return -1;
    y = XInverse(x);
    if (y != NULL)
    rc = XSave(y,oname);
    XFree(x);
    XFree(y);
    Cleanup();

    return rc;
}




/**
@page prog_ziv ziv - Invert

<<<<<<< HEAD
@section syntax Command Line
=======
@section ziv_syntax Command Line
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
<pre>
ziv @em Options @em Input @em Result
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par @em Input
  Input file.
@par @em Result
  Result file.

<<<<<<< HEAD
@section inp Input Files
@par @em Input
  A square natrix or permutation.

@section out Output Files
@par @em Result
  Inverted matrix or permutation.

@section desc Description
=======
@section ziv_inp Input Files
@par @em Input
  A square natrix or permutation.

@section ziv_out Output Files
@par @em Result
  Inverted matrix or permutation.

@section ziv_desc Description
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
This program inverts a matrix or permutation.
**/
