/* ============================= C MeatAxe ==================================
   File:        $Id: zro.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Random orders.
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

static int opt_G = 0;
static int opt_s = 0;
static int Count = 0;
static int NGen = 0;
static void *Gen[MAXGEN];


static MtxApplicationInfo_t AppInfo = { 
"zro", "Random Orders", 
"SYNTAX\n"
"    zro [-GVQs] [-T <#Secs>] <Count> <Gen> ...\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -G ...................... GAP output (implies -Q)\n"
"    -s ...................... Summary only\n"
"\n"
"ARGUMENTS\n"
"    <Count> ................. Number of orders to calculate\n"
"    <Gen> ................... Name of the representation\n"
};


static MtxApplication_t *App = NULL;




/* ------------------------------------------------------------------
   init() - Process command line options and arguments
   ------------------------------------------------------------------ */

static int Init(int argc, const char **argv)

{
    int i;

    App = AppAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;

    /* Process command line options.
       ----------------------------- */
    opt_G = AppGetOption(App,"-G --gap");
    opt_s = AppGetOption(App,"-s --summary");
    if (opt_G) MtxMessageLevel = -100;

    /* Process arguments.
       ------------------ */
    if (AppGetArguments(App,2,MAXGEN+1) < 0)
	return -1;
    Count = atoi(App->ArgV[0]);
    if (Count < 1)
    {
	MTX_ERROR1("Invalid count '%s' (try --help)",App->ArgV[0]);
	return -1;
    }

    /* Read the generators
       ------------------- */
    NGen = App->ArgC - 1;
    for (i = 0; i < NGen; ++i)
    {
	Gen[i] = XLoad(App->ArgV[i+1]);
	if (i > 0 && !XIsCompatible(Gen[0],Gen[i]))
	{
	    MTX_ERROR3("%s and %s: %E",App->ArgV[1],App->ArgV[i+1],
		MTX_ERR_INCOMPAT);
	}
    }

    return 0;
}




static void Cleanup()

{
    AppFree(App);
}



#define MAXORDERS 10
long Order[MAXORDERS];
int CountTab[MAXORDERS];


static void RandomOrders()

{
    int n;
    void *m = NULL;

    m = XDup(Gen[0]);
    MtxRandomInit(0);
    if (opt_G)
	printf("MeatAxe.RandomOrders := [");
    for (n = 0; n < Count; ++n)
    {
	long o;
	o = XOrder(m);
	if (opt_s)
	{
	    int i;
	    for (i = 0; i <= MAXORDERS && CountTab[i] > 0 && Order[i] != o; ++i);
	    if (i < MAXORDERS)
	    {
		Order[i] = o;
		++CountTab[i];
	    }
	}
	else
	{
	    if (n % 15 == 0)
		printf("\n    ");
	    printf("%4ld",o);
    	    if (opt_G && n < Count - 1) 
		printf(",");
	}
	XMul(m,Gen[MtxRandomInt(NGen)]);
    }
    if (opt_G)
	printf("];\n");
    else
    {
	if (opt_s)
	{
	    int i;

	    for (i = 0; i < MAXORDERS && CountTab[i] > 0; ++i) 
	    {
		int k;
		for (k = i + 1; k < MAXORDERS && CountTab[k] > 0; ++k)
		{
		    if (Order[i] > Order[k])
		    {
			long o = Order[i];
			int c = CountTab[i];
			Order[i] = Order[k];
			Order[k] = o;
			CountTab[i] = CountTab[k];
			CountTab[k] = c;
		    }
		}
	    }

	    printf("Order:");
	    for (i = 0; i < MAXORDERS && CountTab[i] > 0; ++i) 
		printf("%6ld",Order[i]);
	    printf("\nCount:");
	    for (i = 0; i < MAXORDERS && CountTab[i] > 0; ++i) 
		printf("%6d",CountTab[i]);
	    printf("\n");
	}
	else
	    printf("\n");
    }
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

    RandomOrders();

    Cleanup();
    return EXIT_OK;
}

/**
@page prog_zro zro - Random Orders

<<<<<<< HEAD
@section syntax Command Line
=======
@section zro_syntax Command Line
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
<pre>
zro [@em Options] [-Gs] @em Count @em Gen [@em Gen ...]
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts.
@par -G
  GAP output.
@par -s
  Summary output.
@par @em Count
  Number or orders to calculate.
@par @em Gen
  Generator.

<<<<<<< HEAD
@section inp Input Files
@par @em Gen
  Generator (matrix or permutation).

@section desc Description
=======
@section zro_inp Input Files
@par @em Gen
  Generator (matrix or permutation).

@section zro_desc Description
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
@b zro calculates the order of @em Count random elements of the
group generated by a set of matrices or permutations. This information
can be helpful to find out which group is generated by a given set of
matrices or permutations.
*/

