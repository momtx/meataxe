/* ============================= C MeatAxe ==================================
   File:        $Id: mkcycl.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     This program calculates a representative for each cyclic
		submodule of the condensed modules.h
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <string.h>
#include <stdlib.h>



/* --------------------------------------------------------------------------
   Global data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO

static MtxApplicationInfo_t AppInfo = { 
"mkcycl", "Find Cyclic Submodules",
"SYNTAX\n"
"    mkcycl [<Options>] <Name>\n"
"\n"
"ARGUMENTS\n"
"    <Name> .................. Name of the representation\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -G ...................... GAP output (implies -Q)\n"
"\n"
"FILES\n"
"    <Name>.cfinfo ........... I Constituent info file\n"
"    <Name><Cf>.{1,2...}k .... I Generators on condensed modules\n"
"    <Name><Cf>.np ........... I Condensed peak words\n"
"    <Name><Cf>.v ............ O Cyclic submodules\n"
};






static MtxApplication_t *App = NULL;


static int opt_G = 0;
static MatRep_t *Rep;		/* The representation */
int NCyclic;			/* Number of cyclic submodules found */
Matrix_t *Cyclic[MAXCYCL];	/* List of cyclic submodules */
Lat_Info LI;			/* Data from .cfinfo file */



/* --------------------------------------------------------------------------
   ParseCommandLine() - Process command line options and arguments
   -------------------------------------------------------------------------- */

static int ParseCommandLine()
{
    opt_G = AppGetOption(App,"-G --gap");
    if (opt_G) 
	MtxMessageLevel = -100;
    if (AppGetArguments(App,1,1) != 1)
	return -1;
    return 0;
}



/* --------------------------------------------------------------------------
   Init() - Program initialization

   This function initializes all global variables.
   -------------------------------------------------------------------------- */

static int Init(int argc, const char **argv)

{
    /* Parse command line
       ------------------ */
    if ((App = AppAlloc(&AppInfo,argc,argv)) == NULL)
	return -1;
    if (ParseCommandLine() != 0)
    {
	MTX_ERROR("Error in command line");
	return -1;
    }

    /* Read the .cfinfo file
       --------------------- */
    MESSAGE(0,("\n*** FIND CYCLIC SUBMODULES ***\n\n"));
    if (Lat_ReadInfo(&LI,App->ArgV[0]) != 0)
    {
	MTX_ERROR1("Error reading %s",App->ArgV[0]);
	return -1;
    }
    return 0;
}



/* --------------------------------------------------------------------------
   Spinup() - Spin up one seed vector

   This function spins up one seed vector and compares the resulting module
   with the list of cyclic submodules found so far. If the module is new, it
   is added to the list.
   -------------------------------------------------------------------------- */

static void Spinup(Matrix_t *seed)

{
    Matrix_t *sub;
    int i;

    sub = SpinUp(seed,Rep,SF_FIRST|SF_SUB,NULL,NULL);
    for (i = 0; i < NCyclic; ++i)
    {
	int issub = IsSubspace(sub,Cyclic[i],1);
	if (issub == -1)
	{
	    MTX_ERROR("Subspace comparison failed");
	    return;
	}
	if (issub && sub->Nor == Cyclic[i]->Nor)
	{
	    MatFree(sub);	/* Module already in list */
	    return;
	}
    }
    if (NCyclic >= MAXCYCL)
    {
	MTX_ERROR1("Too many cyclic submodules (maximum = %d)",MAXCYCL);
	return;
    }
    Cyclic[NCyclic] = sub;
    ++NCyclic;
}


/* --------------------------------------------------------------------------
   WriteResult() - Write the output file

   This function collects the generating vectors of all syclic submodules
   into a single matrix and write this matrix to a file.
   -------------------------------------------------------------------------- */

void WriteResult(int cf, int cond_dim)

{
    int i;
    char fn[100];
    Matrix_t *result;

    result = MatAlloc(FfOrder,NCyclic,cond_dim);
    for (i = 0; i < NCyclic; ++i)
	MatCopyRegion(result,i,0,Cyclic[i],0,0,1,-1);
    sprintf(fn,"%s%s.v",LI.BaseName,Lat_CfName(&LI,cf));
    MESSAGE(1,("Writing %s\n",fn));
    MatSave(result,fn);
    MatFree(result);
}



/* --------------------------------------------------------------------------
   FindCyclic() - Find cyclic submodules

   This function is called once for each constituent. It finds all cyclic
   submodules in the corresponding condensed module, i.e., in the module
   that was condensed with a peak word for this constituent.

   The function uses the seed vector generator to get one vector from each
   one-dimensional subspace. It spins up the seed vectors and stores the 
   distinct submodules obtained in the global variable <Cyclic>.

   Generating vectors for the cyclic sumodules are written to the file
   xxxx.v where xxxx is the constituent name (e.g., test44a).

   Note: 
   Strictly speaking we find all cyclic submodules which are invariant 
   under the condensed generators. The algebra generated by the condensed
   generators may not be the full condensed algebra. For this reason we
   add the kondensed peak word to the set of generators.
   -------------------------------------------------------------------------- */

void FindCyclic(int cf)

{
    Matrix_t *seed;
    Matrix_t *seed_basis;
    char fn[200];
    long vec_no;		/* Seed vector number */
    int count;			/* Number of vectors tried */
    int cond_dim;		/* Dimension of the condensed module */
    int i;

    /* Read the generators and the condensed peak words
       ------------------------------------------------ */
    sprintf(fn,"%s%s.%%dk",LI.BaseName,Lat_CfName(&LI,cf));
    MESSAGE(1,("Loading generators for %s%s\n",LI.BaseName,Lat_CfName(&LI,cf)));
    Rep = MrLoad(fn,LI.NGen);
    sprintf(fn,"%s%s.np",LI.BaseName,Lat_CfName(&LI,cf));
    MrAddGenerator(Rep,MatLoad(fn),0);

    /* Spin up all seed vectors
       ------------------------ */
    cond_dim = Rep->Gen[0]->Nor;
    seed = MatAlloc(FfOrder,1,cond_dim);
    seed_basis = MatId(FfOrder,cond_dim);
    NCyclic = 0;
    count = 0;
    vec_no = 0;
    while ((vec_no = MakeSeedVector(seed_basis,vec_no,seed->Data)) >= 0)
    {
	++count;
	if (count % 100 == 0)
	    MESSAGE(2,("  %d vectors, %d submodules\n",count,NCyclic));
	Spinup(seed);
    }

    /* Write the result and clean up
       ----------------------------- */
    MESSAGE(0,("%s%s: %d cyclic submodule%s (%d vectors tried)\n",
	LI.BaseName,Lat_CfName(&LI,cf),NCyclic,NCyclic == 1 ? " " : "s",
	count));
    WriteResult(cf,cond_dim);
    MatFree(seed);
    MatFree(seed_basis);

    /* Clean up
       -------- */
    for (i = 0; i < NCyclic; ++i)
	MatFree(Cyclic[i]);
    MrFree(Rep);
}



int main(int argc, char *argv[])

{
    int i;

    if (Init(argc,(const char **)argv) != 0)
	return -1;
    for (i = 0; i < LI.NCf; ++i)
	FindCyclic(i);
    if (MtxMessageLevel >= 0)
        printf("\n");
    AppFree(App);
    return 0;
}

/**
@page prog_mkcycl mkcycl - Find Cyclic Subspaces

<<<<<<< HEAD
@section syntax Command Line
=======
@section mkcycl_syntax Command Line
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
<pre>
mkcycl @em Options [-G] @em Name
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -G
  Produce output in GAP format.
@par @em Name
  Name of the representation.

<<<<<<< HEAD
@section inp Input Files
=======
@section mkcycl_inp Input Files
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
@par @em Name.cfinfo
  Constituent info file.
@par @em NameCf.1k, @em NameCf.2k, ...
  Generators on condensed modules.
@par @em NameCf.np
  Condensed peak words

<<<<<<< HEAD
@section out Output Files
@par @em NameCf.v
  Cyclic submodules.

@section desc Description
=======
@section mkcycl_out Output Files
@par @em NameCf.v
  Cyclic submodules.

@section mkcycl_desc Description
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
This program is invoked after @ref prog_pwkond "pwkond" has calculated the
condensation with respect to the peak words. @b mkcycl calculates, for each
condensed module, its 1-dimensional subspaces. The output is a list of vectors
(in matrix form) for each irreducible constituent, which generate all cyclic
submodules. For example, if "X10a" is the constituent's name, the list of
vectors is written to "X10a.v".

<<<<<<< HEAD
@section impl Implementation Details
=======
@section mkcycl_impl Implementation Details
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
@b mkcycl uses a very simple approach: it spins up every vector in the 
condensed module (avoiding scalar multiples, though), and maintains a
list of all cyclic submodules found. As the dimension of the condensed
module grows, the number of vectors to spin up quickly becomes very large. 
This poses an upper limit on the dimension of condensed modules, i.e., on 
the multiplicity of irreducible constituents. Over GF(2), for example, a
16-dimensional condensed module requires about 20 hours of CPU time
on a standard workstation.

A second limit concerns the number of cyclic submodules. Usually there
are much less cyclic submodules than 1-spaces. Sometimes, however, it 
may happen that the peak word found in the second step is "bad" in the 
sense that the condensed generators commute. In such a case one finds 
a large number of cyclic submodules and the following steps will probably 
take too much time. For this reason, the @ref prog_pwkond "pwkond" program
has an option to exclude one or more specified peak words from the search.
So, if the peak word turns out to be "bad", you can try another one.
**/

