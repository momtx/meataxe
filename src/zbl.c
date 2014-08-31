/* ========================== C MeatAxe =============================
   zbl.c - Make matrix lower triangular (keeping bottom left).

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */

#include "meataxe.h"



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static MtxApplicationInfo_t AppInfo = {
    "zbl", "Bottom Left of a Matrix",
    "SYNTAX\n"
    "    zbl " MTX_COMMON_OPTIONS_SYNTAX " <Matrix> <Result>\n"
    "\n"
    "ARGUMENTS\n"
    "    <Matrix> ................ Input matrix\n"
    "    <Result> ................ Output matrix\n"
    "\n"
    "OPTIONS\n"
    MTX_COMMON_OPTIONS_DESCRIPTION
};

static MtxApplication_t *App = NULL;
static const char *iname = 0, *oname = 0;

MTX_DEFINE_FILE_INFO


static int Init(int argc, const char **argv)
{
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;

    /* Parse command line.
       ------------------- */
    if (AppGetArguments(App,2,2) < 0)
	return -1;
    
    iname = App->ArgV[0];
    oname = App->ArgV[1];
    return 0;
}



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char *argv[])
{
    int fl;
    PTR m1;
    int nor, noc, i, j;
    FILE *inp, *out;

    if (Init(argc, (const char **)argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return 1;
    }

    /* Open the input file.
       -------------------- */
    if ((inp = FfReadHeader(iname,&fl,&nor,&noc)) == NULL)
	return 1;
    if (fl < 1)
    {
	MTX_ERROR2("%s: %E",iname,MTX_ERR_NOTMATRIX);
	return 1;
    }

    /* Allocate workspace.
       ------------------- */
    FfSetField(fl);
    FfSetNoc(noc);
    m1 = FfAlloc((long)1);

    /* Open the output file.
       --------------------- */
    if ((out = FfWriteHeader(oname,fl,nor,noc)) == NULL)
    {
	return 1;
    }

    /* Do the job.
       ----------- */
    for (i = 1; i <= nor; ++i)
    {
	FfReadRows(inp,m1,1);
	for (j = i + 1; j <= noc; ++j)
	    FfInsert(m1,j,FF_ZERO);
	FfWriteRows(out,m1,1);
    }

    /* Clean up.
       --------- */
    fclose(inp);
    fclose(out);
    return 0;
}

/**
@page prog_zbl zbl - Bottom Left

<<<<<<< HEAD
@section syntax Command Line
=======
@section zbl_syntax Command Line
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
<pre>
zbl [@em Options] @em Input @em Output
</pre>

@par @em Options
Standard options, see @ref prog_stdopts

@par @em Input
Input matrix

@par @em Output
Result matrix

<<<<<<< HEAD
@section inp Input Files
@par @em Input
Input matrix

@section out Output Files
@par @em Output
Result matrix

@section desc Description
=======
@section zbl_inp Input Files
@par @em Input
Input matrix

@section zbl_out Output Files
@par @em Output
Result matrix

@section zbl_desc Description
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c

This program reads in a matrix, zeroizes all entries above the main diagonal,
and writes out the result. For example:
<pre>
        Mat                  Result
12100201201212212122  10000000000000000000
21000000000000001111  21000000000000000000
12212120101201201201  12200000000000000000
01020010020102012012  01020000000000000000
12100012101201012012  12100000000000000000
12120212121012012012  12120200000000000000
</pre>
Notice that the input matrix need not be square, and the output
matrix has always the same dimensions as the input matrix.

If file names are omitted, the matrix is read from G1
and output goes to P2.

The purpose of this program is to enable the MeatAxe to check if an irreducible representation
in characteristic 2 fixes a quadratic form. This job is not particularly simple -- in many ways
it is just a bodge, but it is possible:
-# Put the representation into standard basis using @ref prog_zsp "zsp" in "standard basis" mode.
-# Find the symplectic form fixed by using @ref prog_zsp "zsp" to make the matrix that conjugates
  the representation to its dual.
-# Quadratic forms can be represented by lower triangular matrices. Since the representation
  is in standard basis (so all the basis vectors are images of the first under the group) the
  diagonal entries of any fixed quadratic form must all be equal, so try each field entry in
  turn by adding that scalar matrix to the bottom half of the symplectic form.
-# For each quadratic form Q made as in 3, test if it is fixed by forming @f$G^TQG@f$ for each
  generator @f$G@f$, and checking that the diagonal is still the same as it was before (the
  symplectic form should have been checked before starting). The check can be done by adding
  the form to the result, then doing @b zbl, @b ztr and @b zbl again --- the result will be the zero
  matrix (use @ref prog_znu "znu") iff the form was fixed (given that the symplectic one was).

*/
