////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Order of a matrix or permutations.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>
#include <string.h>


#define MAXORDER 100000		/* Maximal order */
#define MAXORDER_C 1000		/* Maximal order on cyclic subspaces */


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


static MtxApplicationInfo_t AppInfo = {
"zor", "Order of a Matrix or Permutation",
"SYNTAX\n"
"    zor [-m <MaxOrder>] [-GQVq] <File>\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -m ...................... Set highest possible order\n"
"    -q ...................... Quick mode: Find a lower bound for the order\n"
"    -G ...................... GAP output\n"
"\n"
"FILES\n"
"    <File> .................. I  A matrix or permutation\n"
};


static MtxApplication_t *App = NULL;
static const char *iname;
static MtxFile_t *ifile;
static int maxord = -1;	/* Limit set with -m option */
static int opt_q = 0;		/* Quick mode */
static int opt_G = 0;		/* GAP output */



/* ------------------------------------------------------------------
   ordmat() - Order of a matrix
   ------------------------------------------------------------------ */

static int ordmat()

{   PTR m1, v;
    PTR base, bend;
    int dim;
    int *piv;
    char *ispiv;
    int ord;

    ffSetField(ifile->Field); 
    ffSetNoc(ifile->Noc); 
    if (ifile->Nor != ifile->Noc) 
    {
	mtxAbort(MTX_HERE,"%s: %s",iname,MTX_ERR_NOTSQUARE);
	return 1;
    }
    m1 = ffAlloc(ifile->Noc, ifile->Noc);
    base = ffAlloc(ifile->Noc+1, ifile->Noc);
    piv = NALLOC(int,ifile->Noc);
    ispiv = NALLOC(char,ifile->Noc);
    memset(ispiv,0,ifile->Noc);
    v = ffAlloc(1, ifile->Noc);
    if (mfReadRows(ifile,m1,ifile->Noc) != ifile->Noc)
    {
	mtxAbort(MTX_HERE,"Error reading input file");
	return -1;
    }
    ord = 1;
    bend = base;

    for (dim = 0; dim < ifile->Noc; )
    {	
	PTR start = bend;
	int tord = 0;
	int closed = 0;

	/* Find the next seed vector
	   ------------------------- */
	int i = 0;
	while (i < ifile->Noc && ispiv[i]) ++i;
	MTX_ASSERT(i < ifile->Noc, -1);
	ffMulRow(bend,FF_ZERO);
	ffInsert(bend,i,FF_ONE);

	/* Calculate order on the cyclic subspace
	   -------------------------------------- */
	do
	{   PTR b;
	    FEL f;
	    int pv;

	    /* Save the vector and extend the basis,
	       if the vector is linearly independent.
	       -------------------------------------- */
	    ffCopyRow(v,bend);
	    if (!closed)
	    {	b = base;
	    	for (i = 0; i < dim; ++i)
	    	{   
		    f = ffExtract(bend,piv[i]);
		    if (f != FF_ZERO)
		    {	
			ffAddMulRow(bend,b,ffNeg(ffDiv(f,ffExtract(b,piv[i]))));
		    }
		    ffStepPtr(&b,ifile->Noc);
		}
		pv = ffFindPivot(bend,&f);
		if (pv >= 0)
		{   piv[dim++] = pv;
		    ispiv[pv] = 1;
		    ffStepPtr(&bend, ifile->Noc);
		}
		else
		    closed = 1;
	    }

	    /* Apply the matrix.
	       ----------------- */
	    if (++tord > MAXORDER_C)
	    {  
		mtxAbort(MTX_HERE,"zor: Partial order is over %d",MAXORDER_C);
		return 1;
	    }
	    ffMapRow(v,m1,ffNoc,bend);
	}
	while (ffCmpRows(bend,start));

	/* Calculate l.c.m. of all tord's
	   ------------------------------ */
	for (i = ord; ord % tord != 0; ord += i);
	if (ord > MAXORDER)
	{
	    mtxAbort(MTX_HERE,"zor: Order is over %d",MAXORDER);
	    return 1;
	}
	if (opt_q && dim > ffNoc/10) break;
	if (maxord > 1)
	{   if (ord > maxord)
	    {
	    	fprintf(stderr,"zor: Order is over %d\n",maxord);
	    	exit(1);
	    }
	    if (ord == maxord) break;
	}
    }

    if (opt_q && ord != maxord)
	MESSAGE(0,("ORDER IS A MULTIPLE OF %d\n",ord));
    else
    {
	if (opt_G)
	    printf("MeatAxe.Order := %d;\n",ord);
	else
	    printf("ORDER IS %d\n",ord);
    }
    return 0;
}


/* ------------------------------------------------------------------
   ordperm() - Order of permutations
   ------------------------------------------------------------------ */

static int ordperm()

{
    long order, iper;
    Perm_t *perm;


    if (ifile->Field != -1) 
    {
	mtxAbort(MTX_HERE,"%s: %s",iname,MTX_ERR_NOTPERM);	/* No monomials */
	return -1;
    }
    if ((perm = permAlloc(ifile->Nor)) == NULL)
    {
	mtxAbort(MTX_HERE,"Error allocating permutation");
	return -1;
    }
    if (opt_G) 
	printf("MeatAxe.Orders := [");
    for (iper = 1; iper <= ifile->Noc; ++iper)
    {
	if (mfReadLong(ifile,perm->Data,perm->Degree) != perm->Degree)
    	{
	    mtxAbort(MTX_HERE,"Error reading permutation");
	    return -1;
	}
	Perm_ConvertOld(perm->Data,perm->Degree);
	order = permOrder(perm);
	permFree(perm);
	if (order < 0)
	{
	    mtxAbort(MTX_HERE,"Error calculating order");
	    return -1;
	}
	if (opt_G)
	{
	    if (iper > 1) printf(",");
	    printf("%ld",order);
	}
	else
	    printf("ELEMENT %ld HAS ORDER %ld\n",iper, order);
    }
    if (opt_G) printf("];\n");
    return 0;
}












static int Init(int argc, char **argv)

{
    /* Process command line options.
       ----------------------------- */
    App = appAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;
    opt_G = appGetOption(App,"-G --gap");
    if (opt_G)
	MtxMessageLevel = -100;
    opt_q = appGetOption(App,"-q --quick");
    maxord = appGetIntOption(App,"-m --max-order",-1,1,1000000);

    /* Process arguments.
       ------------------ */
    if (appGetArguments(App,1,1) < 0)
	return -1;
    iname = App->ArgV[0];

    /* Open input file, call the appropriate function
       ---------------------------------------------- */
    if ((ifile = mfOpen(iname)) == NULL)
    {
	mtxAbort(MTX_HERE,"Error opening input file");
	return -1;
    }
    return 0;
}



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char **argv)

{
    int rc;

    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return 1;
    }
    if (ifile->Field < 0)
	rc = ordperm();
    else
	rc = ordmat();
    mfClose(ifile);
    appFree(App);
    return rc;
}


/**
@page prog_zor zor - Order

@section syntax Command Line
<pre>
zor @em Options [-q] [-m @em MaxOrder] @em Input
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -q
  Quick mode (see description).
@par -m @em MaxOrder
  Set upper limit for the oder of cyclic subspaces (see description).
@par @em Input
  Input file.

@section zor_inp Input Files
@par @em Input
  Matrix or permutation.

@section zor_desc Description
This program reads a file, containing either permutations, or a
square matrix, and calculates the order(s) and prints the message
<pre>
ORDER IS xxxx</pre>
There are two options to reduce the run time of the program.
Using the "-m" option you can specify a maximal expected order. 
If, during the algorithm described below, the order reaches 
this limit, the program will stop and print an appropriate message.
The second option, "-q", makes @b zor stop if the dimension of W
(see below) reaches 1/10 of the dimension of the whole space. In 
this case, the message is
<pre>
ORDER IS A MULTIPLE OF @em NNN
</pre>
Note: The "-q" and "-m" options have no effect for permutations.

@section zor_impl Implementation Details
If the input is a matrix, the order is found by calculating the orders on cyclic
subspaces and taking the least common multiple.
The algorithm works as follows:
- Let A be the given matrix and V the space A acts upon.
  Set W:={0} (the trivial subspace) and o:=1.
- (NEXT) Chose a vector v not in W.  Calculate the cyclic subspace
  C generated by v and the order o' of A on C.
- Set o:=lcm(o,o')
- W:=W + C. If W=V, o is the order of A and the program terminates.
  Otherwise, continue with (NEXT).

Gaussian elimination is used to maintain a basis of W in echelon form. 
In order to avoid infinite loops, there is a limit on o'. If the
vector does not return after 1000 multiplications the order is assumed
to be infinite and the program stops with an error message. This happens
also if the value of o exceeds 100000.

If the input file contains permutations, each one is read in and
its order is calculated as the least common multiple of the orbit
sizes. The result is printed in the format
<pre>
ELEMENT nn HAS ORDER nnn
</pre>

The whole matrix plus a second matrix of the same size must fit into
memory. In the case of permutations, there must be enough memory for
one permutation.
**/


// vim:fileencoding=utf8:sw=3:ts=8:et:cin
