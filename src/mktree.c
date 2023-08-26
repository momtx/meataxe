////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Calculate an element tree of a group
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"

#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>





/* --------------------------------------------------------------------------
   Global data
   -------------------------------------------------------------------------- */



/**
 ** A group element.
 **/
typedef struct
{
    Matrix_t *Matrix;
    int Source;
    int Gen;
} Entry_t;

static Entry_t *Elms = NULL;			/* The element list */
static int NElms = 0;				/* Number of elements */
static int MaxNElms = 0;			/* Size of the list */
static const char *Name = NULL;			/* Name of the representation */
static MatRep_t *Rep = NULL;			/* The representation */
static int NoOutput = 0;			/* -n: No output file */

static MtxApplicationInfo_t AppInfo = { 
"mktree", "Enumerate group elements",
"\n"
"SYNTAX\n"
"    mktree [-n] [-g <NGen>] <Name>\n"
"\n"
"ARGUMENTS\n"
"    <Name> .................. Name of the representation\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -g <NGen> ............... Set number of generators (default is 2)\n"
"    -n ...................... Don't write .elt file, print order only\n"
"\n"
"FILES\n"
"    <Name>.{1,2,...} ........ I Generators\n"
"    <Name>.elt .............. O Element tree\n"
};


static MtxApplication_t *App = NULL;




/* --------------------------------------------------------------------------
   IsInList() - Check if a matrix is in our list

   Description:
     This function checks if a given matrix, is already contained in our
	 element list, <Elms>.

   Arguments:
     <mat>: The matrix to be checked.

   Remarks:
     The number of matrices in <Elms> is taken from <TKI.Order>.

   Return:
     1 if the matrix is in <Elms>, 0 otherwise
   -------------------------------------------------------------------------- */

static int IsInList(Matrix_t *mat)
{
    int i;

    for (i = 0; i < NElms; i++)
    {
	if (matCompare(mat,Elms[i].Matrix) == 0)
	    return 1;
    }
    return 0;
}


/* --------------------------------------------------------------------------
   Init() - Initialize

   Description:
     This function initializes all global variables and processes command
	 line options and arguments.

   Arguments:
     <argc>: Number of arguments
	 <argv>: Argument list.
   -------------------------------------------------------------------------- */

static int Init(int argc, char **argv)

{
    int ngen;

    if ((App = appAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;

    /* Command line options.
       --------------------- */
    NoOutput = appGetOption(App,"-n --no-output");
    ngen = appGetIntOption(App,"-g",2,1,1000);

    /* Arguments.
       ---------- */
    if (appGetArguments(App,1,1) < 0)
	return -1;
    Name = App->ArgV[0];

    /* Load the generators.
       -------------------- */
    Rep = mrLoad(Name,ngen);
    if (Rep == NULL)
	return -1;

    /* Initialize the element list.
       ---------------------------- */
    MaxNElms = Rep->NGen + 10;
    if ((Elms = NALLOC(Entry_t,MaxNElms)) == NULL)
    {
	mtxAbort(MTX_HERE,"Cannot allocate lement list");
	return -1;
    }
    NElms = 0;


    return 0;
}




/* --------------------------------------------------------------------------
   AddToList() - Add one element to the element tree

   Description:
     This function adds one element to the element tree.

   Arguments:
	 <mat>: The new element.
	 <src>: Index of element <mat> was calculated from-
	 <gen>: Number of generator (1,2,3,...) used to calculate <mat> from
		<elem[src]>.

   Remarks:
     Both <src> and <gen> start counting with 1, not 0!
   -------------------------------------------------------------------------- */

static int AddToList(Matrix_t *mat, int src, int gen)

{
    /* Extend the list, if necessary.
       ------------------------------ */
    if (NElms >= MaxNElms)
    {
	int newmax = MaxNElms + 100;
	Entry_t *newlist = NREALLOC(Elms,Entry_t,newmax);
	if (newlist == NULL)
	{
	    mtxAbort(MTX_HERE,"Cannot extend element list, req size=%d",newmax);
	    return -1;
	}
	Elms = newlist;
	MaxNElms = newmax;
    }
    Elms[NElms].Matrix = mat;
    Elms[NElms].Source = src;
    Elms[NElms].Gen = gen;
    ++NElms;
    return 0;
}


static int MakeTree()

{
    int rc = 0;
    int src;

    if (AddToList(matId(ffOrder,ffNoc),-1,-1) != 0)
	return -1;

    /* Calculate all elements
       ---------------------- */
    for (src = 0; src < NElms; ++src)
    {
	int g;
	for (g = 0; g < Rep->NGen; ++g)
	{
	    /* Calculate next element
	       ---------------------- */
	    Matrix_t *newelem = matDup(Elms[src].Matrix);
	    matMul(newelem,Rep->Gen[g]);

	    /* If it is new, add to tree, else discard
	       --------------------------------------- */
	    if (IsInList(newelem))
	    {
		matFree(newelem);
		if (src == 0)
		    MESSAGE(0,("Warning: generator %d is redundant\n",g+1));
	    }
	    else
	    {
	    	MESSAGE(2,("%d x %d = %d\n",src,g,NElms));
		if (AddToList(newelem,src,g) != 0)
		{
		    rc = -1;
		    break;
		}
		if (NElms % 50 == 0)
		    MESSAGE(0,("%d elements\n",NElms));
	    }
	}
    }
    MESSAGE(0,("Done. The group has %d elements.\n",NElms));
    return rc;
}

static void WriteOutput()

{
    char fn[200];
    IntMatrix_t *mat;
    int i;

    sprintf(fn,"%s.elt",Name);
    MESSAGE(1,("Writing %s\n",fn));
    mat = imatAlloc(NElms,2);
    if (mat == NULL)
	return;
    for (i = 0; i < NElms; ++i)
    {
	mat->Data[2*i] = Elms[i].Source;
	mat->Data[2*i + 1] = Elms[i].Gen;
    }
    imatSave(mat,fn);
    imatFree(mat);
}



static void Cleanup()

{
    appFree(App);
}

/* --------------------------------------------------------------------------
   main() - Program entry point
   -------------------------------------------------------------------------- */


int main(int argc, char **argv)

{
    int rc = 0;

    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return 1;
    }

    if (MakeTree() != 0)
	rc = 1;
    if (rc == 0 && !NoOutput)
	WriteOutput();

    Cleanup();
    return rc;
}




/**
@page prog_mktree mktree - Enumerate Group Elements

@section mktree_syntax Command Line
<pre>
mktree @em Options [-n] [-g @em NGen] @em Name
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -n
  Do not write the outpt file.
@par -g
  Set the number of generators (default: 2).
@par @em Name
  Name of the representation.

@section mktree_inp Input Files
@par @em Name.1, @em Name.2, ...
  Generators.

@section mktree_out Output Files
@par @em Name.elt
  Element tree.

@section mktree_desc Description
This program enumerates all elements of a finitely generated matrix
group. By default, the program assumes that the group has two
generators, which are read from @em Name.1 and @em Name.2.
A different number of generators can be specified with "-g".

Unless the "-n" option is used, the program writes the element tree
to @em Name.elt. The element tree describes how the group elements
can be calculated as products of generators. It is actually
a matrix with two columns and one row for each group element.
The i-th row of this matrix describes how the i-th element is 
calculated:

- The row (-1,-1) represents the unit element. This row appears
  first in the output.
- (0,k) means that this element is the k-th generator. Note that
  the generator number starts with 0, i.e., the first generator
  has k=0.
- (s,k) means that the corresponding element is obtained by
  multiplying the s-th element from the right by the
  k-th generator. Both generator and element numbers start
  with 0. Thus, the (0,k) lines decribed above are actually
  a special case of the general (s,k) lines.

The following example shows the output for a group of the order 10
with two generators, a and b:

<pre>
Line   Contents     Meaning
----   --------     ---------------
 1      -1   -1     identity matrix
 2      0    0      a
 3      0    1      b
 4      1    0      aa
 5      1    1      ab
 6      3    0      aaa
 7      3    1      aab
 8      5    0      aaaa
 9      5    1      aaab
10      7    1      aaaab
</pre>

@section mktree_impl Implementatino Details
The program holds all group elements in memory. This limits the application
of the program to fairly small groups and representations of small degree.
**/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
