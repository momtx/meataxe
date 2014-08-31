/* ============================= C MeatAxe ==================================
   File:        $Id: zvp.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Vector permute (make permutations from matrices)
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include <stdlib.h>
#include <string.h>
#include "meataxe.h"



#define MAXVEC 100000	/* Default value for maxvec */

#define MAX_GENERATORS	50


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

static MtxApplicationInfo_t AppInfo = { 
"zvp", "Vector Permute",
"SYNTAX\n"
"    zvp [<Options>] [-g <NGen>] <Mat> <Seed> <Perm> [<Orbit>]\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -g <NGen> ............... Set number of generators (default: 2)\n"
"    -n ...................... No output\n"
"    -p ...................... Permute 1-spaces instead of vectors\n"
"    -v ...................... Write Vectors to <Orbit>\n"
"    -l <Limit> .............. Set maximal orbit size\n"
"    -m ...................... Make (generate) seed vectors from <Seed>\n"
"    -s <N> .................. Start with seed vector <N>\n"
"\n"
"ARGUMENTS\n"
"    <Mat> ................... Generator base name\n"
"    <Seed> .................. Seed vector file name\n"
"    <Perm> .................. Output file name\n"
"    <Orbit> ................. Orbit file name, default os 'orbit'\n"
"\n"
"FILES\n"
"    <Mat>.{1,2...} .......... I Generators (square matrices)\n"
"    <Seed> .................. I Seed vectors (matrix)\n"
"    <Perm>.{1,2...} ......... O Permutations\n"
"    <Orbit> ................. O The orbit (matrix)\n"
};
static MtxApplication_t *App = NULL;


static int NGen = 2;			/* Number of generators */
static Matrix_t *Gen[MAX_GENERATORS];	/* Generators */
static Matrix_t *Seed;			/* Seed space */
static int Generate = 0;		/* Generate vectors (-m) */
static int SeedStart = 0;		/* First seed vector to use (-s) */
static int maxvec = MAXVEC;		/* Max. orbit size (-l) */

static int tabsize;			/* Size of hash table */
static int *vecpos;		        /* Position of vectors in table */
static int *vecno;			/* Numbers of vectors in table */
static char *isfree;
static int nvec;		        /* Number of vectors obtained so far */
static int nfinished;			/* Number of finished vectors */
static PTR vtable;			/* Holds the vectors */
static PTR tmp;
static long *perm[MAX_GENERATORS];	/* Permutations */


static long iseed = 0;			/* Current seed vector */
int proj = 0;				/* Operate on 1-spaces (-p)  */
static int vecout = 0;			/* Output vectors, too (-v) */
static int noout = 0;			/* No output, print orbit sizes only (-n) */
static long hasm, hasp, neh;		/* Parameters for the hash function */

static const char *MatName;		/* Matrix name */
static const char *PermName;		/* Permutation name */
static const char *SeedName;		/* Seed space name */
static const char *OrbName = "orbit";	/* Orbit file name (optional) */





/* ------------------------------------------------------------------
   InitHash() - Initialize the hash parameters
   ------------------------------------------------------------------ */

static void InitHash()

{	
    int dim = Gen[0]->Noc;
    int tod, i;

    neh = (dim > 25) ? 25 : dim;
    tod = 100;
    hasp = tabsize;
    if (hasp < 110) 
	tod = 7;
    if (hasp > 11)
    {	
	i = 1;
	while (++i <= tod)
	{	
	    if (hasp % i == 0)
	    {	
		--hasp;
		i = 1;
	    }
	}
    }
    if (hasp < 100)
	hasm = 3;
    else
    {	
	if (FfChar % 83 == 0)
	    hasm = 89;
	else
	    hasm = 83;
    }
}



static int ParseCommandLine()

{
    /* Command line options.
       --------------------- */
    noout = AppGetOption(App,"-n");
    vecout = AppGetOption(App,"-v");
    proj = AppGetOption(App,"-p");
    NGen = AppGetIntOption(App,"-g",2,1,MAX_GENERATORS);
    maxvec = AppGetIntOption(App,"-l",MAXVEC,0,-1);
    Generate = AppGetOption(App,"-m");
    SeedStart = AppGetIntOption(App,"-s",1,1,10000000) - 1;

    /* Command line arguments.
       ----------------------- */
    if (AppGetArguments(App,3,4) < 0)
	return -1;
    MatName = App->ArgV[0];
    SeedName = App->ArgV[1];
    PermName = App->ArgV[2];
    if (App->ArgC > 3)
	OrbName = App->ArgV[3];
    return 0;
}


static int ReadFiles()

{
    int i;

    /* Read the generators.
       -------------------- */
    for (i = 0; i < NGen; ++i)
    {
	char fn[200];
	sprintf(fn,"%s.%d",MatName,i+1);
	if ((Gen[i] = MatLoad(fn)) == NULL)
	    return -1;
	if (Gen[i]->Nor != Gen[i]->Noc)
	{
	    MTX_ERROR2("%s: %E",fn,MTX_ERR_NOTSQUARE);
	    return -1;
	}
	if (i > 0 && (Gen[i]->Field != Gen[0]->Field 
	    || Gen[i]->Nor != Gen[0]->Nor))
	{
	    MTX_ERROR4("%s and %s.%d: %E",fn,MatName,i+1,MTX_ERR_INCOMPAT);
	    return -1;
	}
    }

    /* Read the seed space.
       -------------------- */
    if ((Seed = MatLoad(SeedName)) == NULL)
	return -1;
    if (Seed->Field != Gen[0]->Field || Seed->Noc != Gen[0]->Nor)
    {
	MTX_ERROR3("%s and %s.0: %E",SeedName,MatName,MTX_ERR_INCOMPAT);
	return -1;
    }

    return 0;
}


static int AllocateTables()

{
    int i;

    tabsize = maxvec + maxvec / 10 + 1;
    MESSAGE(1,("Allocating tables (size=%d)\n",tabsize));
    vtable = FfAlloc(tabsize+1);
    tmp = FfAlloc(1);
    vecpos = NALLOC(int,tabsize+1);
    vecno = NALLOC(int,tabsize+1);
    isfree = NALLOC(char,tabsize+1);
    for (i = 0; i < NGen; ++i)
    {
	perm[i] = NALLOC(long,maxvec+1);
	if (perm[i] == NULL) 
	{
	    MTX_ERROR("Error allocationg permutation");
	    return -1;
	}
    }

    InitHash();
    return 0;
}


static int Init(int argc, const char **argv)

{   
    App = AppAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;
    if (ParseCommandLine() != 0)
	return -1;
    if (ReadFiles() != 0)
	return -1;
    if (AllocateTables() != 0)
	return -1;
    iseed = SeedStart;
    return 0;
}





/* ------------------------------------------------------------------
   Normalize() - Normalize a vector
   ------------------------------------------------------------------ */

static void Normalize(PTR row)

{	
    FEL f;

    FfFindPivot(row,&f);
    FfMulRow(row,FfInv(f));
}



/* ------------------------------------------------------------------
   MakeNextSeedVector() - Make the next seed vector in <tmp>
   ------------------------------------------------------------------ */

static int MakeNextSeedVector()

{
    MESSAGE(1,("Starting with seed vector %ld\n",iseed));
    if (Generate)
    {
	iseed = MakeSeedVector(Seed,iseed,tmp);
	if (iseed < 0)
	    return -1;
    }
    else
    {
	if ((int) iseed >= Seed->Nor)
	    return -1;
	FfCopyRow(tmp,MatGetPtr(Seed,(int)iseed));
	++iseed;
    }
    if (proj)
	Normalize(tmp);
    return 0;
}


/* ------------------------------------------------------------------
   hash() - The hash function
   ------------------------------------------------------------------ */

static int hash(PTR row)

{	
    register int pos = 0, i;

    for (i = 0; i < neh; ++i)
	pos = (pos * hasm + FfToInt(FfExtract(row,i))) % hasp;
    return pos;
}





/* ------------------------------------------------------------------
   InitTables() - Prepare everything for spin-up. Assumes that the
   seed vector is in <tmp>.
   ------------------------------------------------------------------ */

static void InitTables()

{
    int i;
    int pos;
    PTR x;

    for (i = 0; i < (int) tabsize; ++i)
      isfree[i] = 1;
    nvec = 1;
    nfinished = 0;
    pos = hash(tmp);
    x = FfGetPtr(vtable,pos);
    FfCopyRow(x,tmp);
    isfree[pos] = 0;
    vecpos[0] = pos;
    vecno[pos] = 0;
}


/* ------------------------------------------------------------------
   mkorbit() - Make the orbit
   ------------------------------------------------------------------ */

static int MakeOrbit()

{
    int pos, im, pos1;
    PTR x;
    int igen;	/* Apply which generator */

    igen = 0;
    while (nfinished < nvec && nvec <= maxvec)
    {
	MESSAGE(3,("Vec[%d] * Gen[%d] = ",nfinished,igen));
	x = FfGetPtr(vtable,vecpos[nfinished]);
	FfMapRow(x,Gen[igen]->Data,FfNoc,tmp);
	if (proj) 
	    Normalize(tmp);

	/* Look up the vector in the hash table.
	   ------------------------------------- */
	pos1 = pos = hash(tmp);
	x = FfGetPtr(vtable,pos);
	while (!isfree[pos] && FfCmpRows(tmp,x))
	{
	    if (++pos == tabsize)
	    {	
		pos = 0;
		x = vtable;
	    }
	    else
	    	FfStepPtr(&x);
	    /* The following should never happen since the
	       hash table is always larger than maxvec */
	    MTX_VERIFY(pos != pos1);
	}

	if (isfree[pos])	/* New vector */
	{
	    MESSAGE(3,("%d (new)\n",nvec));
	    FfCopyRow(x,tmp);
	    im = nvec;
	    isfree[pos] = 0;
	    vecpos[nvec] = pos;
	    vecno[pos] = nvec;
	    nvec++;
	}
	else
	{
	    im = vecno[pos];
	    MESSAGE(3,("%d\n",im));
	}
	if (nvec % 10000 == 0)	
	    MESSAGE(2,("%d vectors, %d finished\n",nvec,nfinished));
	perm[igen][nfinished] = im;

	/* Next generator
	   -------------- */
	if (++igen >= NGen)
	{
	    igen = 0;
	    ++nfinished;
	}
    }

    if (nfinished < nvec)
	return -1;
    return 0;
}




static void Cleanup()

{
    AppFree(App);
}




static int WriteOutput()

{
    int rc = 0;
    int i;
    char fn[200];

    if (noout) 
	return 0;

    /* Write vectors
       ------------- */
    if (vecout)
    {	
	MtxFile_t *f;
	int i;

	MESSAGE(1,("Writing orbit to %s\n",OrbName));
	if ((f = MfCreate(OrbName,FfOrder,nvec,FfNoc)) == 0)
	{
	    MTX_ERROR1("Cannot open %s",OrbName);
	    rc = -1;
	}
        for (i = 0; i < nvec; ++i)
	    {	
		PTR x = FfGetPtr(vtable,vecpos[i]);
		MfWriteRows(f,x,1);
	    }
        MfClose(f);
    }

    /* Write permutations
       ------------------ */
    MESSAGE(1,("Writing permutations\n"));
    for (i = 0; i < NGen; ++i)
    {   
	MtxFile_t *f;
	sprintf(fn,"%s.%d",PermName,i+1);
	/* if ((f = MfCreate(OrbName,-1,nvec,1)) != NULL) */
	if ((f = MfCreate(fn,-1,nvec,1)) == NULL)
	{
	    MTX_ERROR1("Cannot open %s",fn);
	    rc = -1;
	    continue;
	}
	if (MfWriteLong(f,perm[i],nvec) != nvec)
	{
	    MTX_ERROR("Cannot write data");
	    rc = -1;
	    continue;
	}
	MfClose(f);
    }
    return rc;
}




/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, const char **argv)

{   
    int rc = 1;

    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return 1;
    }


    /* Main loop
       --------- */
    while (1)
    {
	if (MakeNextSeedVector() != 0)
	{
	    MTX_ERROR("No nore seed vectors");
	    break;
	}
	InitTables();
	if (MakeOrbit() == 0)
	{
	    MESSAGE(0,("Vector %ld: Orbit size = %d\n",iseed,nvec));
	    WriteOutput();
	    rc = 0;
	    break;
	}
        MESSAGE(0,("Orbit of vector %ld is longer than %d\n",iseed,maxvec));
    }
    Cleanup();
    return rc;
}



/**
@page prog_zvp zvp - Vector permute

<<<<<<< HEAD
@section syntax Command Line
=======
@section zvp_syntax Command Line
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
<pre>
zvp [@em Options] [-g @em NGen] [-l @em Limit] [-s @em Start] [-np] @em Man @em Seed @em Perm [@em Orbit]
</pre>

@par @em Options
Standard options, see @ref prog_stdopts.

@par -g @em NGen
    Set the number generators (default is 2).

@par -n
    No output to files, messages only.

@par -p
    "Projective" mode, permute 1-spaces instead of vectors.
    
@par -v
    Write the orbit to @em Orbit.

@par -l @em Limit
    Set the orbit size limit. Default is 100000.

@par -m
    Generate seed vectors by takig all possible linear combinations of the rows of
    @em Seed.  Without this option, only the rows of @em Seed are used as seed vectors.

@par -s @em Start
    Start with the given seed vector number instead of 1.

@par @em Mat
    Name of the representation.

@par @em Seed
    Seed space file name.

@par @em Perm
    Permutation file name.

@par @em Orbit
    Orbit file name. Default is "orbit".

<<<<<<< HEAD
@section inp Input Files
=======
@section zvp_inp Input Files
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c

@par @em Mat.1, @em Mat.2, ...
    Generators (square matrices).

@par @em Seed
    Seed vectors (matrix).

<<<<<<< HEAD
@section out Output Files
=======
@section zvp_out Output Files
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c

@par @em Perm.1, @em Perm.2, ...
    Permutations.

@par @em Orbit
    The orbit (matrix).

<<<<<<< HEAD
@section desc Description
=======
@section zvp_desc Description
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
This program reads a set of matrices and one or more vectors from 
@em Seed, and finds the orbit of the vector under the matrices. The 
action of the matrices on the orbit is written out in permutation form.

By default, two matrices are read from @em Mat.1 and @em Mat.2.
You can specify a different number of matrices using the "-g" option,
but the naming scheme is always the same. All matrices must be square, 
over the same field, and of equal dimension. The seed space, @em Seed, 
must be a matrix over the same field and the number of columns must match the matrices.

The program tries seed vectors until no more seed vectors are
available or the orbit is small enough (as specified by "-l").
With "-p" all vectors in the orbit file are normalized, i.e., their first
non-zero entry is equal to one.

<<<<<<< HEAD
@subsection impl Implementation Details
=======
@subsection zvp_impl Implementation Details
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
After initializing everything and reading the input files the program
enters the main loop. The next seed vector is read in, or generated,
and the orbit is set up initially to contain only this vector.
Then, the complete orbit is calculated by applying the 
matrices to vectors in the orbit until no new vectors appear,
or until the limit is reached. 
In the latter case, the program proceeds with the next seed vector.
If an orbit has been found, the action of the two matrices on the 
orbit is written out in permutation format an appropriate message
is printed. If the `-v' option was used, also the orbit is written 
out.

The matrices, seed vectors, and all vectors in the orbit must fit
into memory.

**/
