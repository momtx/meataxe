////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Make orbits under permutations
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <stdlib.h>
#include <string.h>


#define MAXPERM 50		/* Max. number of permutations */
#define STACKSIZE 100000


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


static Perm_t *Perm[MAXPERM];			/* Permutations */
static const char *orbname;
static const char *permname;
static int nperm = 2;
static int Degree;
static int Seed = 0;
static int32_t *OrbNo;
static uint32_t *OrbSize;
static int NOrbits;
static long Stack[STACKSIZE];
static int Sp = -1;


static MtxApplicationInfo_t AppInfo = { 
"zmo", "Make Orbits",
"SYNTAX\n"
"    zmo [<Options>] [-g <#Perms>] <Perm> <Orbits>\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -g <#Perms> ............. Set number of permutations (default: 2)\n"
"    -s <Seed> ............... Set seed point (default: 1)\n"
"\n"
"FILES\n"
"    <Perm>.{1,2...} ......... I Permutations\n"
"    <Orbits> ................ O Orbit table and sizes\n"
};
static MtxApplication_t *App = NULL;



////////////////////////////////////////////////////////////////////////////////////////////////////

static void readPermutations()
{
    for (int i = 0; i < nperm; ++i)
	Perm[i] = permLoad(strEprintf("%s.%d",permname,i+1));
    Degree = Perm[0]->degree;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    nperm = appGetIntOption(App,"-g",2,1,MAXPERM);
    Seed = appGetIntOption(App,"-s --seed",1,1,1000000) - 1;
    appGetArguments(App,2,2);
    permname = App->argV[0];
    orbname = App->argV[1];
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void cleanup()
{
   for (int i = 0; i < nperm; ++i)
      permFree(Perm[i]);
    if (OrbNo!= NULL) 
	sysFree(OrbNo);
    appFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void allocWorkspace()
{
    OrbNo = NALLOC(int32_t,Degree);
    for (int i = 0; i < Degree; ++i)
	OrbNo[i] = -1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int makeOrbits()
{
    int points_remaining;
    int seedpos = 0;

    MTX_LOGD("Finding orbits");
    Stack[++Sp] = Seed;
    OrbNo[Seed] = 0;
    NOrbits = 1;
    for (points_remaining = Degree; points_remaining > 0; --points_remaining)
    {
	long pt;
	long orb;
	int i;

	// If there is somthing on the stack, take it. Otherwise start a new orbit.
	if (Sp >= 0)
	{
	    pt = Stack[Sp--];
	    orb = OrbNo[pt];
	}
	else
	{
	    for (pt = seedpos; pt < Degree && OrbNo[pt] >= 0; ++pt)
		;
	    if (pt >= Degree)
	    {
		mtxAbort(MTX_HERE,"Internal error: no point found to continue");
		return -1;
	    }
	    seedpos = pt + 1;
	    orb = NOrbits++;
	    OrbNo[pt] = orb;
	}

	// Apply all permutations.
	for (i = 0; i < nperm; ++i)
	{
	    long image = Perm[i]->data[pt];
	    if (OrbNo[image] < 0)
	    {
		OrbNo[image] = orb;
		if (Sp >= STACKSIZE - 1)
		{
		    mtxAbort(MTX_HERE,"Stack overflow");
		    return -1;
		}
		Stack[++Sp] = image;
	    }
	    else
	    {
		if (OrbNo[image] != orb)
		{
		    mtxAbort(MTX_HERE,"Internal error: inconsistent orbit numbers");
	    	    return -1;
		}
	    }
	}
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void CalcSizes()
{
    MTX_LOGD("Calculating orbit sizes");
    OrbSize = NALLOC(uint32_t,NOrbits);
    memset(OrbSize,0,sizeof(*OrbSize) * NOrbits);
    for (int i = 0; i < Degree; ++i)
	++OrbSize[OrbNo[i]];
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void WriteOutput() 
{
    FILE *f = sysFopen(orbname,"wb");
    // Write the orbit number table.
    {
    uint32_t hdr[3] = {MTX_TYPE_INTMATRIX, 1, Degree};
    sysWrite32(f,hdr,3);
    sysWrite32(f, OrbNo, Degree);
    }

    // Write the orbit sizes table.
    {
    uint32_t hdr[3] = {MTX_TYPE_INTMATRIX, 1, NOrbits};
    sysWrite32(f,hdr,3);
    sysWrite32(f,OrbSize,NOrbits);
    }
    fclose(f);

    // Print orbit sizes
    int size[20];
    int count[20];
    int dis = 0;
    int i;
    for (i = 0; i < NOrbits; ++i)
    {
       int k;
       for (k = 0; k < dis && size[k] != OrbSize[i]; ++k);
       if (k < dis)
          ++count[k];
       else if (dis < 20)
       {
          size[dis] = OrbSize[i];
          count[dis] = 1;
          ++dis;
       }
    }
    for (i = 0; i < dis; ++i)
       MTX_LOGI("%6d ORBIT%c OF SIZE %6d", count[i],count[i] > 1 ? 'S':' ',size[i]);
}



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char **argv)

{
    init(argc,argv);
    readPermutations();
    allocWorkspace();
    makeOrbits();
    CalcSizes();
    WriteOutput();
    cleanup();
    return (EXIT_OK);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_zmo zmo - Make Orbits

@section zmo_syntax Command Line
<pre>
zmo @em Options [-g @em NPerms] [-s @em Seed] @em Perm @em Orbits
</pre>

@par @em Options
Standard options, see @ref prog_stdopts

@par [-g @em NPerms]
Set the number of permutatins (default is 2).

@par [-s @em Seed]
Set the seed point (default is 1).

@par @em Perm
Permutation base name.

@par @em Orbits
Output file name.

@section zmo_inp Input Files
@par @em Perm.1, @em Perm.2, ...
Permutations.

@section zmo_out Output Files
@par @em Orbits
Orbit number and sizes.


@section zmo_desc Description

This program calculates the orbits under a set of permutations.
By default, the program works with two permutations which are read
from @em Perm.1 and @em Perm.2.  Using the `-g' option you
can specify a different number of permutations (using the same file 
name convention). For example,
<pre>
zmo -g 3 p orbs
</pre>
reads three permutations from p.1, p.2, and p.3.
All permutations must have the same degree.

The result is written to @em Orbits and consists of two parts:
- The orbit number table. This is an integer matrix with one
  row and N columns, where N is the degree of the permutations.
  The i-th entry in this table contains the orbit number of
  point i. Note: Orbit numbers start with 0.
- The orbit sizes table. This is an integer matrix with one
  row and K columns, where K is the number of orbits. The 
  i-th entry contains the size of the orbit number i.
This file can be fed into @ref prog_zkd "zkd" or @ref prog_zuk "zuk".
At the end, the program prints a message containing the orbit sizes.
Note that at most 20 different orbit sizes are shown here.


@section zmo_impl Implementation Details
The algorithm uses a fixed size stack to store points. At the beginning,
a seed point is searched which has not yes assigned an orbit number.
This point is assigned the next orbit number (beginning with 0) and put on the stack.
The main part consists of taking a point from the stack and applying all generators.
Any new points obtained in this way are assigned the same orbit 
number and put on the stack. This step is repeated until the stack is
empty, i.e., the orbit is exhausted. Then, the next seed point is
seached, and its orbit is calculated, and so on, until all orbits
are found.
By default, the first seed point is 1. A different seed point may
be selected with the "-s" option.

The number of permutations must be less than 50. All permutations must 
fit into memory at the same time. Also the stack size is limited to 
100000 positions.
*/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
