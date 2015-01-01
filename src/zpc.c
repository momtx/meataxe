////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Permutation chop.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <string.h>
#include <stdlib.h>
#include "meataxe.h"


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

static MtxApplicationInfo_t AppInfo = {
"zpc", "Permutation Chop",
"SYNTAX\n"
"    zpc [-b] <Perm1> <Perm2> <Seed> <S1> <S2> <Q1> <Q2>\n"
"    zpc [-b] -g <#Perm>[.<#Gen>] <Perm> <Seed> <Sub> <Quot>\n"
"\n"
"OPTIONS\n"
"    -b ...................... Find action on block system\n"
"    -g ...................... Multiple permutation mode\n"
"\n"
"FILES\n"
"    <Perm1,2> ............... I Generators: permutations of equal degree\n"
"    <Seed> .................. I Seed: One point or (with -b) one block\n"
"    <S1,2> .................. O Action on the orbit or block system\n"
"    <Q1,2> .................. O Action on the rermaining points\n"
"\n"
};



static MtxApplication_t *App = NULL;
long *seed;			/* Read from G3 */
int npoints;			/* Degree of input permutations */
int blksize;			/* Block size */
int nblocks;			/* Number of blocks */
long *perm[MAXGEN+1];		/* Generators */
int ngen = 2;			/* Number of generators */
int nperm = 2;			/* Number of permutations to chop */
long *orb;
long *num;
long start;
long orblen;

const char *genname[MAXGEN+1];
const char *subname[MAXGEN+1];
const char *quotname[MAXGEN+1];
const char *seedname;

int opt_b = 0;			/* Option -b (permute blocks) */
int opt_g = 0;			/* Option -g (additional generators) */


static char *helptext[] = {
NULL};





/* ------------------------------------------------------------------
   mkname()
   ------------------------------------------------------------------ */

static char *mkname(const char *basename,int i)

{   
    char *c;

    c = NALLOC(char,strlen(basename)+10);
    sprintf(c,"%s%d",basename,i);
    return c;
}


/* ------------------------------------------------------------------
   setnperm()
   ------------------------------------------------------------------ */

static void setnperm(const char *c)

{
    nperm = -1;
    opt_g = 1;
    if (sscanf(c,"%d.%d",&nperm,&ngen) == 1)
	ngen = nperm;
    if (nperm > MAXGEN || ngen > nperm || ngen < 2)
	MTX_ERROR1("-g: %E",MTX_ERR_OPTION);
    if (ngen != nperm && MSG1)
        MTX_ERROR2("%d generators, %d permutations",ngen,nperm);
}


/* ------------------------------------------------------------------
   parseargs() - Process command line options and arguments
   ------------------------------------------------------------------ */

static int parseargs()

{
    int i;
    const char *c;
    int args_needed;

    opt_b = AppGetOption(App,"-b --block-system");
    c = AppGetTextOption(App,"-g",NULL);

    opt_g = c != NULL;
    if (opt_g)
	setnperm(c);

    /* Get file names
       -------------- */
    args_needed = opt_g ? 4 : 7;
    if (AppGetArguments(App,args_needed,args_needed) < 0)
	return -1;
    if (opt_g)
    {
	seedname = App->ArgV[1];
	for (i = 1; i <= nperm; ++i)
	{	
	    genname[i] = mkname(App->ArgV[0],i);
	    subname[i] = mkname(App->ArgV[2],i);
	    quotname[i] = mkname(App->ArgV[3],i);
	}
    }
    else
    {	
	quotname[2] = App->ArgV[6];
	quotname[1] = App->ArgV[5];
	subname[2] = App->ArgV[4];
	subname[1] = App->ArgV[3];
	seedname = App->ArgV[2];
	genname[2] = App->ArgV[1];
	genname[1] = App->ArgV[0];
    }
    return 0;
}



/* ------------------------------------------------------------------
   init() - Initialize everything, read input files
   ------------------------------------------------------------------ */

static int init(int argc, char **argv)

{
    long *m;
    int i;
    MtxFile_t *f;


    App = AppAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;
    parseargs(argc, argv);

    /* Read the Permutation
       -------------------- */
    for (i = 1; i <= nperm; ++i)
    {
	if ((f = MfOpen(genname[i])) == NULL)
	    return -1;
	if (f->Field != -1 || f->Noc != 1)
	{
	    MTX_ERROR2("%s: %E",genname[i],MTX_ERR_NOTPERM);
	    return -1;
	}
	if (i == 1)
	    npoints = f->Nor;
	else if (f->Nor != npoints)
	{
	    MTX_ERROR3("%s and %s: %E",genname[1],genname[i],MTX_ERR_INCOMPAT);
	    return -1;
	}
    	m = NALLOC(long,npoints);
	if (m == NULL)
	    return -1;
	if (MfReadLong(f,m,npoints) != npoints)
	    return -1;
	Perm_ConvertOld(m,npoints);
        perm[i] = m;
	MfClose(f);
    }

    /* Seed
       ---- */
    if ((f = MfOpen(seedname)) == NULL)
	return -1;
    if (f->Field != -1) 
    {
	MTX_ERROR3("%s: %E (found type %d, expected -1)",seedname,
	    MTX_ERR_FILEFMT,(int)f->Field);
    }
    if (!opt_b)
    {	
	blksize = 1;
	nblocks = npoints;
    }
    else
    {	
	blksize = f->Nor;
	nblocks = npoints / blksize;
	if (npoints % blksize != 0)
	{
	    MTX_ERROR("BLOCK SIZE DOES NOT DIVIDE DEGREE");
	    return -1;
	}
    }

    seed = NALLOC(long,blksize);
    if (seed == NULL) 
	return -1;
    if (MfReadLong(f,seed,blksize) != blksize)
	return -1;
    MfClose(f);

    /* Allocate tables
       --------------- */
    if ((orb = NALLOC(long,npoints+1)) == NULL 
	|| (num = NALLOC(long,npoints+1)) == NULL
       )
       return -1;
    for (i = 0; i < npoints; ++i)
    {
	orb[i] = -1;
	num[i] = -1;
    }

    return 0;
}


/* ------------------------------------------------------------------
   writeresult() - Calculate action on orbits and write output files.
   ------------------------------------------------------------------ */

static int writeresult()

{   
    int cosize;
    long *s, *p;
    int i;
    long k, im;
    MtxFile_t *f;

    cosize = npoints - orblen;
    if (cosize == 0)
	printf("Transitive on %ld %s\n",orblen/blksize,
		opt_b ? "blocks" : "points");
    else
	printf("Intransitive - 'sub' %ld  'quot' %ld\n",
		orblen/blksize,cosize/blksize);

    /* Calculate action on first orbit
       ------------------------------- */
    if (orblen % blksize != 0) 
    {
	MTX_ERROR2("Invalid block system: orblen=%d, blksize=%d",
	    orblen,blksize);
    }
    s = NALLOC(long,orblen/blksize+1);
    for (i = 1; i <= nperm; ++i)
    {	
	p = perm[i];
	for (k = 0; k < orblen/blksize; ++k)
	    s[k] = -1;
	for (k = 0; k < npoints; ++k)
	{
	    long x, y;
	    if (num[k] >= orblen) continue;
	    x = num[k] / blksize;
	    y = num[p[k]] / blksize;
	    if (s[x] == -1)
		s[x] = y;
	    else
	    {
		if (s[x] != y) 
		    MTX_ERROR("Invalid block system");
	    }
	}
	if ((f = MfCreate(subname[i],-1,orblen/blksize,1)) == NULL)
	    return -1;
	if (MfWriteLong(f,s,orblen/blksize) != (orblen/blksize))
	    return -1;
	MfClose(f);
    }
    free(s);
    if (opt_b || cosize <= 0) 
	return -1;

    /* Calculate action on other orbits (not with -b option)
       ----------------------------------------------------- */
    if ((s = NALLOC(long,cosize+1)) == NULL)
	return -1;
    for (i = 1; i <= nperm; ++i)
    {	
	p = perm[i];
	for (k = 1; k < npoints; ++k)
	{
	    if (num[k] < orblen) 
		continue;
	    im = p[k];
	    s[num[k]-orblen] = num[im]-orblen;
	}
	if ((f = MfCreate(quotname[i],(long)-1,cosize,(long)1)) == NULL)
	    return -1;
	if (MfWriteLong(f,s,cosize) != cosize)
	    return -1;
	MfClose(f);
    }
    free(s);
    return 0;
}


/* ------------------------------------------------------------------
   chop() - Make orbit
   ------------------------------------------------------------------ */

static void chop()

{	
    long level, i, k, newpt;
    int gen;

    /* Set up the orbit to contain only the seed point/block
       ----------------------------------------------------- */
    for (i = 0; i < blksize; ++i)
    {	
	orb[i] = seed[i];
	num[orb[i]] = i;
    }
    orblen = blksize;
    level = 0;

    /* Make orbit
       ---------- */
    while (level < orblen)
    {	/* Apply all generators
           -------------------- */
	for (gen = 1; gen <= ngen; ++gen)
	{	
	    if (num[perm[gen][orb[level]]] == -1)
			/* New point ? */
	    {	
		for (i = 0; i < blksize; ++i)
		{	
		    newpt=perm[gen][orb[level+i]];
		    num[newpt] = orblen+i;
		    orb[orblen+i] = newpt;
		}
		orblen += blksize;
	    }
	}
	level += blksize;
    }

    if (opt_b)	/* No `quotient' when permuting blocks */
	return;

    /* Renumber the remaining points
       ----------------------------- */
    k = orblen;
    for (i = 0; i < npoints; ++i)
    {	
	if (num[i] == -1)
	    num[i] = ++k;
    }
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char *argv[])

{
    if (init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return 1;
    }

    chop();
    writeresult();
    AppFree(App);
    return EXIT_OK;
}


/**
@page prog_zpc zpc - Permutation Chop

@section zpc_syntax Command Line
<pre>
zpc [@em Options] [-b] @em Perm1 @em Perm2 @em Seed @em S1 @em S2 @em Q1 @em Q2

zpc [@em Options] [-b] -g @em NPerm[.@em NGen] @em Perm @em Seed @em Ssub @em Quot
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -b
  Find action on block system.
@par -g @em -g @em NPerm[.@em NGen]
  Set the number of Permutations.
@par @em Perm1, @em Perm2
  Generators (without -g)
@par @em Perm
  Generator base name (with -g)
@par @em Seed
  Seed point or block.
@par @em S1, @em S2
  Action on the orbit or block system (without -g).
@par @em Sub
  Base name for the action on the orbit or block system (with -g).
@par @em Q1, @em Q2
  Action on the remaining points (without -g).
@par @em Quot
  Base name for the action on the remaining points (with -g).
  

@section zpc_inp Input Files
@par @em Perm1, @em Perm2
  Generators (without -g).
@par @em Perm.1, @em Perm.2, ...
  Generators (with -g).
@par @em Seed
  Seed point or block.

@section zpc_out Output Files
@par @em S1, @em S2
  Action on the orbit or block system (without -g).
@par @em Sub.1, @em Sub.2, ...
  Action on the orbit or block system (with -g).
@par @em Q1, @em Q2
  Action on the remaining points (without -g).
@par @em Quot.1, @em Quot.2, ...
  Action on the remaining points (with -g).

@section zpc_desc Description
If invoked without any options, this program reads two permutations from @em Perm1 and @em Perm2
and a point from @em Seed. The orbit containing that point is made, and the action on the orbit
is written to @em S1 and @em S2.
If any points remain, the action on these is written to @em Q1 and @em Q2.

With the "-b" option, @em  Seed is taken as a list of points defining one block of a
block system for the permutation group, and the program computes the action of
both the generators on the block system.
Nothing is written to @em Q1 and @em Q2 in this case,
even if the group generated by the two permutations does not act
transitively on the block system.

A number of permutations different from 2 can be specified by using the
"-g" option. @em NPerm is the number of permutations to chop, and @em NGen
specifies how many of these are sufficient to generate the whole permutation group.
For example, with "-g 5.2" @b zpc would chop 5 permutations, assuming that the first
two generate the whole group.
@em NGen may be omitted, the default being @em NPerm.
Thus, instead of "-g 5.5" you may simply type "-g 5".
When the "-g" option is used, there must be exactly 4 Arguments, 3
of which (@em Perm, @em Sub and @em Quot) are interpreted as "base names".
For example,
<pre>
zpc -g 3  z seed s q
</pre>
would read three permutations from "z1", "z2" and "z3",
and a point from "seed". The actions on the orbit would be
written to "s1", "s2" and "s3", and the actions on the
remaining points to "q1", "q2" and "q3".
The number of permutations, specified with "-g", must not be greater than 10.

The program checks the input files to see if they contain reasonable data.
Failure of any of these checks produces an error message and the program stops.
If all goes well, the first point of the permutation in @em Seed is
taken to be the first point in the orbit. The orbit is then calculated
by applying the generators to all new points.
The points of the orbit are renumbered, and, for each of the
permutations, the action on the orbit is written out. If this orbit is
all the points, the message
<pre>
TRANSITIVE ON {\em n} POINTS</pre>
appears, and the program stops. If any points remain, these are also
renumbered, and the action on the remaining points is written out, too
(unless "-b" was used). In this case
<pre>
INTRANSITIVE - 'SUB' n  'QUOT' m
</pre>
is written to stdout, the 'SUB' being the orbit of the given point.

With the "-b" option, the orbit is calculated in the same way, but
a whole block, i.e., a set of points, is read from @em Seed and
taken as the starting "point". The points are then renumbered in such
a way that the orbit containing the seed is {1,2,...,N, where
N is a multiple of the block size, and the blocks are contiguous.
After this rearrangement, point number n belongs to block ⌊(n-1)/b⌋,
where b is the block size. The action of the permutations on the blocks
can then be calculated by simply dividing the point numbers by the block size.
Of course, the block size must be a divisor of the permutation's degree.
Otherwise the program stops with
<pre>
ZPC ERROR - BLOCK SIZE DOES NOT DIVIDE DEGREE
</pre>

It might turn out during the computation, that two points in one
block are mapped to images in different blocks. This means that
@em Seed was not a block of any block system for the given
representation and the program will stop with the error message
<pre>
ZPC ERROR - NOT A BLOCK SYSTEM
</pre>
If all goes well, the user is informed with one of the messages
<pre>
TRANSITIVE ON n BLOCKS
INTRANSITIVE - 'SUB' n  'QUOT' m
</pre>
and the result is written to @em S1 and @em S2.
@em Q1 and @em Q2 are never written when "-b" is specified,
even if the action on blocks is intransitive.

The use if this program is twofold. Firstly it can be used to obtain
transitive constituents from intransitive permutation representations.
Secondly, if the starting point is in some sense canonical, the
resulting renumbered representation will be according to a `standard'
numbering of the points. It should be noticed that this comment only
applies to the representation on the 'SUB' — that written to @em S1 and @em S2.
If this standard-base-like property of the output is not required,
there is no particular reason to choose @em Seed, and point 1,
or indeed one of the permutations can be used.

**/

