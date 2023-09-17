////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Make a permutation which maps A to B
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


#define MAXPERM 50		/* Max. number of permutations */
static Perm_t *Perm[MAXPERM];			/* Permutations */
static long Degree;
static int nperm;
static int Seed;
static int Stop;
static const char *permname;
static const char *scriptname;
static long *ptnr;
static long *pre;
static int *gen;

static MtxApplicationInfo_t AppInfo = { 
"orbrep", "Find word mapping seed to stop point",
"SYNTAX\n"
"    orbrep [<Options>] [-g <#Perms>] <Perm> <Seed> <Stop> <Script>\n"
"\n"
"ARGUMENTS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -g <#Perms> ............. Set number of permutations (default: 2)\n"
"    <Seed> .................. Seed point (1..N)\n"
"    <Stop> .................. Stop point (1..N)\n"
"\n"
"FILES\n"
"    <Perm>.{1,2...} ......... I Permutations\n"
"    <Script> ................ O Word\n"
};
static MtxApplication_t *App = NULL;


/* ###################################################################### */
static int ReadPermutations()
{
    int i;
    char fn[200];
    for (i = 0; i < nperm; ++i)
    {
	sprintf(fn,"%s.%d",permname,i+1);
	Perm[i] = permLoad(fn);
	if (Perm[i] == NULL)
	    return -1;
	if (i > 1 && Perm[i]->Degree != Perm[0]->Degree)
	{
	    mtxAbort(MTX_HERE,"%s and %s.1 have different degrees",
		fn,permname);
	    return -1;
	}
    }
    Degree = Perm[0]->Degree;
    if (Seed < 0 || Seed >= Degree)
    {
	mtxAbort(MTX_HERE,"Illegal seed point, valid range is 1..%ld.",Degree+1);
	return -1;
    }
    if (Stop < 0 || Stop >= Degree)
    {
	mtxAbort(MTX_HERE,"Illegal stop point, valid range is 1..%ld.",Degree+1);
	return -1;
    }

    return 0;
}

/* ###################################################################### */
static int Init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;

    /* Command line.
       ------------- */
    nperm = appGetIntOption(App,"-g",2,1,MAXPERM);
    if (appGetArguments(App,4,4) < 0)
	return -1;
    permname = App->ArgV[0];
    scriptname = App->ArgV[1];
    Seed = atoi(App->ArgV[2]) - 1;
    Stop = atoi(App->ArgV[3]) - 1;
    return 0;
}

/* ###################################################################### */
static void Cleanup()
{
  if (ptnr != NULL)
    {
      sysFree(ptnr);
      sysFree(pre);
      sysFree(gen);
    }
    appFree(App);
}

/* ###################################################################### */
static int AllocWorkspace()
{
    int i;
    ptnr = NALLOC(long,Degree);
    pre  = NALLOC(long,Degree);
    gen  = NALLOC(int,Degree);
    if (ptnr == NULL || pre == NULL || gen == NULL)
      return -1;
    for (i = 0; i < Degree; ++i)
      {
	ptnr[i] = -1;
	pre[i] = -1;
	gen[i] = -1;
      }
    return 0;
}

/* ###################################################################### */
static int MakeOrbit()
{
    long orblen = 0;
    long done = -1;
    long pt, i, image;

    MESSAGE(1,("Finding orbit of seed point %d\n",Seed+1));
    if (Seed == Stop)
    {
	mtxAbort(MTX_HERE,"Stop point equals seed point");
	return -1;
    }
    ptnr[orblen] = Seed;
    pre[Seed] = Seed;
    while (done < orblen && orblen <= Degree-1)
    {
        pt = ptnr[++done];
	for (i = 0; i < nperm; ++i)
	{
	    image = Perm[i]->Data[pt];
 	    if (pre[image] < 0)
	    {
		ptnr[++orblen] = image;
		pre[image] = pt;
		gen[image] = i;
	    }
	    if (image == Stop)
	    {
		MESSAGE(1,("Stop point %d found\n",Stop+1));
		return 0;
	    }
	}
    }
    mtxAbort(MTX_HERE,"Stop point %d not in orbit",Stop+1);
    return -1;
}

/* ###################################################################### */
static int WriteOutput() 

{
    FILE *f;
    int i = 0;
    int len = 0;
    long pt;
    long *found;

    pt = Stop;
    while (pre[pt] != pt)
      {
	pt = pre[pt];
        len++;
      }
    found = NALLOC(long,len);
    pt = Stop;
    while (pre[pt] != pt)
      {
        found[i++] = pt;
	pt = pre[pt];
      }
    if ((f = sysFopen(scriptname,"w")) == NULL)
	return -1;
    fprintf(f,"word:=[\n");
    for (i = len-1; i >= 0; --i)
    {
      fprintf(f,"%d",gen[found[i]]+1);
      if (i>0)
	fprintf(f,",");
    }
    fprintf(f,"];\n");
    fclose(f);
    sysFree(found);
    return 0;
}

/* ###################################################################### */
int main(int argc, char **argv)
{
    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return 1;
    }
    if (ReadPermutations() != 0)
    {
	mtxAbort(MTX_HERE,"Error reading input files");
	return 1;
    }
    if (AllocWorkspace() != 0)
    {
	mtxAbort(MTX_HERE,"Error allocating workspace");
	return 1;
    }
    if (MakeOrbit() != 0)
    {
	mtxAbort(MTX_HERE,"Error making orbit");
	return 1;
    }
    if (WriteOutput() != 0)
	return 1;
    Cleanup();
    return (EXIT_OK);
}


/**
@page prog_orbrep orbrep - Find a permutation which maps A to B

@section orbrep_syntax Command Line
<pre>
orbrep [@em Options] [-g @em NPerms] @em Perm @em Seed @em Stop @em Script
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -g @em NPerms
  Number of permutations (default is 2).
@par @em Perm
  Permutation base name.
@par @em Seed
  Seed point to start the orbit algorithm.
@par @em Stop
  Stop point for the orbit algorithm.
@par @em Script
Result (see description).

@section orbrep_inp Input Files
@par @em Perm.1, @em Perm.2, ...
  Permutations.

@section orbrep_out Output Files
@par @em Script
  Result.

@section orbrep_desc Description
Given a set of generating permutations, this program find a product of the
generators which maps a given point, @em Start, to a second given
point, @em Stop.
By default, the program works with two generators which are read
from @em Perm.1 and @em Perm.2, respectively.
You can specify a different number of permutations using the -g option.

The result is a list of numbers defining the result as a product
in the generators. This file is a text file an can be read by GAP.

@see @ref prog_zvp

*/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
