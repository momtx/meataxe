////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - This program calculates the null space of a matrix.
////////////////////////////////////////////////////////////////////////////////////////////////////



#include "meataxe.h"
#include <stdlib.h>



/* ------------------------------------------------------------------
   Variables
   ------------------------------------------------------------------ */


static MtxApplicationInfo_t AppInfo = { 
"znu", "Matrix Null-Space", 
"SYNTAX\n"
"    znu [-GQVn] <Matrix> [<NullSpace>]\n"
"\n"
"ARGUMENTS\n"
"    <Matrix> ................ Input file name\n"
"    <Nullspace> ............. Output file name\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -G, --gap ............... GAP output\n"
"    -n, --no-echelon ........ Do no convert the null-space to echelon form\n"
"\n"
};

static MtxApplication_t *App = NULL;
static int opt_G = 0;
static int opt_n = 0;				/* -n: no echelon form */
static const char *matname = NULL;
static const char *nspname = NULL;
static Matrix_t *Matrix = NULL;



static int Init(int argc, char **argv)
{
    /* Process command line options.
       ----------------------------- */
    App = appAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;
    opt_G = appGetOption(App,"-G --gap");
    opt_n = appGetOption(App,"-n --no-echelon");
    if (opt_G)
	MtxMessageLevel = -100;

    /* Process arguments.
       ------------------ */
    if (appGetArguments(App,1,2) < 0)
	return -1;
    matname = App->argV[0];
    if (App->argC > 1)
	nspname = App->argV[1];

    /* Read the matrix.
       ---------------- */
    MESSAGE(1,("Reading %s\n",matname));
    Matrix = matLoad(matname);
    if (Matrix == NULL)
	return 1;

    return 0;
}



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char **argv)

{   
    int nspdim;

    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return 1;
    }

    if (nspname != NULL)
    {
	Matrix_t *null_space = matNullSpace_(Matrix,opt_n ? 1 : 0);
	if (null_space == NULL)
	    return 1;
        MESSAGE(1,("Writing null-space to %s\n",nspname));
	matSave(null_space,nspname);
	nspdim = null_space->nor;
    }
    else
    {
	int old_nor = Matrix->nor;
	matEchelonize(Matrix);
	nspdim = old_nor - Matrix->nor;
    }


    if (opt_G)
    {
	printf("MeatAxe.Nullity := %d;\n",nspdim);
	fflush(stdout);
    }
    else
        MESSAGE(0,("NULLITY %d\n",nspdim));


    appFree(App);
    return 0;
}





/**
@page prog_znu znu - Null-Space

@section znu_syntax Command Line
<pre>
znu @em Options [-G] @e Matrix [@em NullSpace]
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -G
  Output in GAP format.
@par @em Matrix
  The matrix.
@par @em NullSpace
  The null-space.
  
@section znu_inp Input Files
@par @em Matrix
  The matrix (M×N).

@section znu_out Output Files
@par @em NullSpace
  The null-space, a L×M matrix in echelon form.

@section znu_desc Description
This program reads in a matrix and outputs a basis for its null-space in
echelon form. If the @em Nullspace argument is omitted the null-space is
not written out, but its dimension is still printed.

Notice that the input matrix does not need to be square.

@section znu_impl Implementation Details
After reading the matrix, the program generates the n×n identity matrix in
memory where n is the number of rows. It then proceeds to perform row operations
on the matrix until it is in echelon form. The same row operations are
performed on the identity matrix, and whenever a row in the original matrix
becomes zero, the corresponding row of the other matrix is marked for output.
The null-space is always reduced to echelon form.
*/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
