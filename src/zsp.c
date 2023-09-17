////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Spinup, split, and standard basis
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <string.h>
#include <stdlib.h>



/* --------------------------------------------------------------------------
   Global data
   -------------------------------------------------------------------------- */


static Matrix_t *Seed = NULL;		/* Seed vectors */
static MatRep_t *Rep = NULL;		/* Generators */
static Matrix_t *Span = NULL;		/* Invariant subspace */
static IntMatrix_t *OpTable = NULL;	/* Spin-up script */
static int Dim = -1;                    // Dimension / permutation size
static int ngen = 2;			/* Number of generators */
static int opt_G = 0;		        /* GAP output */
static long SeedVecNo = 0;		/* Current seed vector */
static const char *GenName[MAXGEN];	/* File name for generators */
static const char *SeedName = NULL;	/* File name for seed vectors */
static const char *SubspaceName = NULL;	/* File name for invariant subspace */
static const char *SubName = NULL;	/* File name for action on subspace */
static const char *QuotName = NULL;	/* File name for action on quotient */
static const char *OpName = NULL;	/* File name for spin-up script */
static int Permutations = 0;		/* Spin up with permutations */
static const Perm_t *Perm[MAXGEN];	/* Permutations */
static int MaxDim = -1;			/* -d: subspace dim. limit */
static int TryOneVector = 0;		/* -1: try only one seed vector */
static int TryLinearCombinations = 0;	/* -m: try all linear combinations */
static int FindCyclicVector = 0;	/* -e: find a cyclic vector */
static int FindClosure = 0;		/* -c: make closure of the seed space */
static int Standard = 0;		/* -t: standard basis */
static int MaxTries = 0;		/* -x: max. # of tries */



static MtxApplicationInfo_t AppInfo = { 
"zsp", "Spinup, split, and standard basis",
"\n"
"SYNTAX\n"
"    zsp [<Options>] <Gen1> <Gen2> <Seed>\n"
"    zsp [<Options>] [-g <#Gen>] <Gen> <Seed>\n"
"\n"
"ARGUMENTS\n"
"    <Gen1>, <Gen2> .......... Generator names\n"
"    <Gen> ................... Generator name (with -g)\n"
"    <Seed> .................. Seed vector file name\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -b <Basis> .............. Output a basis of the invariant subspace\n"
"    -s <Sub> ................ Calculate the action on the subspace\n"
"    -q <Quot> ............... Calculate the action on the quotient\n"
"    -o <Script> ............. Write a spin-up script\n"
"    -G ...................... GAP output (implies -Q)\n"
"    -g <#Gen> ............... Set number of generators\n"
"    -n <Num> ................ Start with vector <Num>\n"
"    -d <Dim> ................ Set an upper limit for the subspace dimension\n"
"    -1 ...................... Try only one seed vector\n"
"    -m ...................... Make (generate) seed vectors\n"
"    -e ...................... Find a cyclic vector\n"
"    -c ...................... Combine, make the closure\n"
"    -t ...................... Make standard basis\n"
"    -x <Max> ................ Assume subspace is closed after <Max>\n"
"                              multiplications without finding new vector.\n"
"\n"
"FILES\n"
"    <Gen1>,<Gen2>............ I  Generators (without -g)\n"
"    <Gen>.{1,2...} .......... I  Generators (with -g)\n"
"    <Seed> .................. I  Seed vectors\n"
"    <Sub>.{1,2...} .......... O  Action on the subspace (with -s)\n"
"    <Quot>.{1,2...} ......... O  Action on the quotient (with -q)\n"
"    <Basis> ................. O  Basis of the invariant subspace (with -b)\n"
"    <Script> ................ O  Spin-up script (with -o)\n"
};

static MtxApplication_t *App = NULL;


/* ------------------------------------------------------------------
   ReadGenerators() - Read the generators into <Rep>
   ------------------------------------------------------------------ */

static void readGenerators()
{
   int i;

   {
      MtxFile_t *f = mfOpen(GenName[0]);
      mfReadHeader(f);
      uint32_t objectType = mfObjectType(f);
      if (objectType == MTX_TYPE_PERMUTATION)
         Permutations = 1;
      if (objectType == MTX_TYPE_MATRIX) {
         Permutations = 0;
         ffSetField(f->header[0]);
      }
      else {
         mtxAbort(MTX_HERE, "%s: unsupported object type 0x%lx",
               GenName[0], (unsigned long) objectType);
      }
      mfClose(f);
   }

   if (Permutations)
   {
      if (SubName != NULL || QuotName != NULL)
         mtxAbort(MTX_HERE,"'-s' and '-q' are not supported for permutations");
      for (i = 0; i < ngen; ++i)
      {
         Perm[i] = permLoad(GenName[i]);
      }
      Dim = Perm[0]->Degree;
   }
   else
   {
      Rep = mrAlloc(0,NULL,0);
      for (i = 0; i < ngen; ++i)
      {
         Matrix_t *gen = matLoad(GenName[i]);
         mrAddGenerator(Rep,gen,0);
      }
      Dim = Rep->Gen[0]->Noc;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void readSeed()
{
    MtxFile_t *sf = mfOpen(SeedName);
    mfReadHeader(sf);
    if (mfObjectType(sf) != MTX_TYPE_MATRIX)
	mtxAbort(MTX_HERE,"%s: %s",SeedName,MTX_ERR_NOTMATRIX);
    const uint32_t field = sf->header[0];
    const int norSeed = sf->header[1];
    const int nocSeed = sf->header[2];

    if (Permutations)
	ffSetField(field);
    if (nocSeed != Dim || sf->header[0] != ffOrder)
	mtxAbort(MTX_HERE,"%s and %s: %s",GenName[0],SeedName,MTX_ERR_INCOMPAT);
    
    size_t skip = 0;
    if (!TryLinearCombinations && SeedVecNo > 0) {
       // skip the first «SeedVecNo» rows
       if ((skip = SeedVecNo - 1) > norSeed)
          skip = norSeed;
       sysFseekRelative(sf->File, skip * ffRowSize(nocSeed));
    }
    uint32_t num_seed = TryOneVector ? 1 : norSeed - skip;

    Seed = matAlloc(ffOrder,num_seed,Dim);
    mfReadRows(sf,Seed->Data,norSeed,Dim);
    mfClose(sf);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init_args()
{
    int i;

    SubspaceName = appGetTextOption(App,"-b --basis",NULL);
    SubName = appGetTextOption(App,"-s --subspace-action",NULL);
    QuotName = appGetTextOption(App,"-q --quotient-action",NULL);
    OpName = appGetTextOption(App,"-o --script",NULL);

    MaxDim = appGetIntOption(App,"-d --dimension-limit",-1,1,100000000);
    MaxTries = appGetIntOption(App,"-x --max-tries",-1,1,1000);
    TryOneVector = appGetOption(App,"-1 --single-shot");
    TryLinearCombinations = appGetOption(App,"-m --seed-generate");
    FindCyclicVector = appGetOption(App,"-e --find-cyclic-vector");
    FindClosure = appGetOption(App,"-c --combine");
    Standard = appGetOption(App,"-t --standard-basis");

    opt_G = appGetOption(App,"-G --gap");
    ngen = appGetIntOption(App,"-g",-1,1,MAXGEN);
    SeedVecNo = appGetIntOption(App,"-n",0,1,10000000);

    if (MaxDim != -1 && (FindClosure || FindCyclicVector))
	mtxAbort(MTX_HERE,"'-d' cannot be used together with '-c' or '-e'");

    if (   FindClosure
	&& (FindCyclicVector || TryOneVector || TryLinearCombinations)
	)
    if (FindCyclicVector && FindClosure)
    {
	mtxAbort(MTX_HERE,"'-c' cannot be combined with any of '-e', '-1', '-m'");
    }


    if (TryLinearCombinations && TryOneVector)
    {
	mtxAbort(MTX_HERE,"Options '-n' and '-1' cannot be used together");
    }
    if (TryLinearCombinations && SeedVecNo > 0)
    {
	mtxAbort(MTX_HERE,"Options '-m' and '-n' cannot be used together");
    }


    if (opt_G)
	MtxMessageLevel = -100;

    /* Process arguments (generator and seed names).
       --------------------------------------------- */
    if (ngen == -1)
    {
	ngen = 2;
	appGetArguments(App,3,3);
	GenName[0] = App->ArgV[0];
	GenName[1] = App->ArgV[1];
	SeedName =  App->ArgV[2];
    }
    else
    {
	char buf[200];
	appGetArguments(App,2,2);
	SeedName =  App->ArgV[1];
	for (i = 0; i < ngen; ++i)
	{
	    char *c;
	    sprintf(buf,"%s.%d",App->ArgV[0],i+1);
	    GenName[i] = c = sysMalloc(strlen(buf)+1);
	    strcpy(c,buf);
	}
    }
}



static int WriteAction()

{
    MatRep_t *sub = NULL, *quot = NULL;
    MatRep_t **subp = SubName != NULL ? &sub : NULL;
    MatRep_t **quotp = QuotName != NULL ? &quot : NULL;

    if (Split(Span,Rep,subp,quotp) != 0)
    {
	mtxAbort(MTX_HERE,"Split failed");
	return -1;
    }
    if (SubName != NULL)
    {
	mrSave(sub,SubName);
	mrFree(sub);
    }
    if (QuotName != NULL)
    {
	mrSave(quot,QuotName);
	mrFree(quot);
    }
    return 0;
}


static int WriteResult()
{
    if (Span->Nor < Dim && (Standard || FindCyclicVector))
	MESSAGE(0,("ZSP: Warning: Span is only %d of %d\n",Span->Nor,Dim));
    else if (Span->Nor == Dim && TryLinearCombinations)
	MESSAGE(0,("ZSP: Warning: No invariant subspace found\n"));
    else
    {
	MESSAGE(0,("Subspace %d, quotient %d\n",Span->Nor, Span->Noc - Span->Nor));
    }

    /* Write the invariant subspace.
       ----------------------------- */
    if (SubspaceName != NULL)
	matSave(Span,SubspaceName);
    
    /* Write <Op> file
       --------------- */
    if (OpName != NULL)
	imatSave(OpTable,OpName);

    /* Write the action on the subspace and/or quotient.
       ------------------------------------------------- */
    if (SubName != NULL || QuotName != NULL)
	WriteAction();

    return 0;
}




static void init(int argc, char **argv)
{
    App = appAlloc(&AppInfo,argc,argv);
    init_args();
    readGenerators();
    readSeed();
}


static void Cleanup()

{
    appFree(App);
}



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */


int main(int argc, char **argv)
{
    IntMatrix_t **optab = NULL;
    int flags;
    SpinUpInfo_t SpInfo;

    init(argc,argv);

    if (OpName != NULL)
	optab = &OpTable;

    SpinUpInfoInit(&SpInfo);
    if (MaxDim > 0)
	SpInfo.MaxSubspaceDimension = MaxDim;
    if (MaxTries > 0)
	SpInfo.MaxTries = MaxTries;

    flags = 0;

    /* Seed mode: SF_MAKE, SF_EACH, or SF_FIRST.
       ----------------------------------------- */
    if (TryOneVector)
	flags |= SF_FIRST;
    else if (TryLinearCombinations)
	flags |= SF_MAKE;
    else
	flags |= SF_EACH;

    /* Search mode: SF_CYCLIC, SF_SUB, or SF_COMBINE.
       ---------------------------------------------- */
    if (FindCyclicVector)
	flags |= SF_CYCLIC;
    else if (FindClosure)
	flags |= SF_COMBINE;
    else
	flags |= SF_SUB;

    /* Spin-up mode: SF_STD or nothing.
       -------------------------------- */
    if (Standard)
	flags |= SF_STD;

    if (!Permutations)
	Span = SpinUp(Seed,Rep,flags,optab,&SpInfo);
    else
	Span = SpinUpWithPermutations(Seed,ngen,Perm,flags,optab,&SpInfo);
    if (Span == NULL)
    {
	mtxAbort(MTX_HERE,"Spin-up failed");
	return -1;
    }
    if (WriteResult() != 0)
	return 1;
    Cleanup();
    return 0;
}


/**
@page prog_zsp zsp - Spin Up

@section zsp_syntax Command Line
<pre>
zsp [@em Options] [-1emc] [-b @em Bas] [-s @em Sub] [-q @em Quot] [-o @em Scr] [-n @em Vector]
    [-d @em MaxDim] [-x @em MaxTries] @em Gen1 @em Gen2 @em Seed

zsp [@em Options] [-1emc] [-b @em Bas] [-s @em Sub] [-q @em Quot] [-o @em Scr] [-n @em Vector]
    [-d @em MaxDim] [-x @em MaxTries] [-g @em NGen] @em Gen @em Seed
</pre>

@par @em Options
Standard options, see @ref prog_stdopts
@par -b @em Bas
  Write a basis of the invariant subspace to @em Bas.
@par -c
  Combine the span of all seed vectors.
@par -s @em Sub
  Calculate the action on the subspace and write the matrices to @em Sub.1, @em Sub.2, ...
@par -q @em Quot
  Calculate the action on the quotient and write the matrices to @em Quot.1, @em Quot.2, ...
@par -o @em Scr
    Write a spin-up script to @em Scr.
@par -d @em MaxDim
    Set an upper limit of the subspace dimension.
@par -n @em Vector
    Start with the specified seed vector.
@par -1
    Spin up only one vector.
@par -e
    Find a cyclic vector.
@par -p
    Find a proper subspace.
@par -m
    Generate seed vectors.
@par -t
    Spin up canonically (standard basis).
@par -x
    Assume the subspace is closed, if @em MaxTries vectors have been multiplied by
    all generators without yielding a new vector.
@par -g 
    Set the number of generators.
@par -G
    GAP output.
@par @em Gen1
  First generator, if "-g" is not used.
@par @em Gen2
  Second generator, if "-g" is not used.
@par @em Gen
  Base name for generators with "-g".
@par @em Seed
  Seed space file name.


@section zsp_inp Input Files
@par @em Gen1
  First generator, if "-g" is not used.
@par @em Gen2
  Second generator, if "-g" is not used.
@par @em Gen.1, @em Gen.2, ...
  Generators (with "-g").
@par @em Seed
  Seed vectors.

@section zsp_out Output Files
@par @em Sub.1, @em Sub.2, ...
  Action on the subspace (with -s).
@par @em Quot.1, @em Quot.2, ...
  Action on the quotient (with -q).
@par @em Basis.
  Basis of the invariant subspace (with -b).
@par @em Script
  Spin-up script (with -o).

@section zsp_desc Description
This program takes as input a set of matrices or permutations (the "generators"), 
and a list of seed vectors. It uses the spin-up algorithm to find a subspace which
is invariant under the generators. If the generators are matrices, @b zsp can optionally split the representation, i.e., calculate the
action of the generators on both subspace and quotient. Splitting is currently
not possible for permutations.


@subsection inp Specifying Input Files
There are two ways to inkove @b zsp. The first form, without `-g' expects three 
arguments, the two generators and the seed vector file. For example,
||    zsp mat1 mat2 seed
reads the generators from `mat1' and `mat2', and the seed vector from 
`seed'.

If the number of generators is not two, you must use the second form, 
which expects only two arguments. The first argument is treated as a base name.
The actual file names are built by appending suffixes ".1", ".2",... to @em Gen.
For example,
<pre>
zsp -g 3 module seed</pre>
reads three genrators from "module.1", "module.2", and "module.3".
The last argument, @em Seed is always treated as a single file name, 
containing the seed vectors. Of course, the seed vectors must be 
compatible with the generators, i.e., they must be over the save field
and have the same number of columns. The generators must be square 
matrices of the same size and over the same field.

@subsection seedmode Specifying the Seed Mode
@b zsp has three ways of interpreting the seed vector file. The default is
to treat @em Seed as a list of seed vectors, which are used ony-by-one
until one seed vector is successful (see below for the meaning of 
successful), or until all vectors have been used. Normally, @b zsp starts
with the first row of @em Seed, but this can be changed using the `-n' 
option. For example,
<pre>
zsp -n 4 gen1 gen2 seed</pre>
starts with spinning up the fourth row of "seed". If this is not 
successful, @b zsp continues with row 5 and so on up to the end of the
seed vector file.

With "-1" @b zsp spins up only the first seed vector and stops, even if 
the spin-up was not successful. You can use "-n" to select a different
row as seed vector. If any of these options is used, @b zsp loads only
the seed vectors that are actually needed.

If you use the "-m" option, @b zsp treats @em Seed as the basis of a 
seed space and tries all 1-dimensional subspaces as seed vectors. In
this mode, seed vectors are constructed by taking linear combinations
of the rows of @em Seed. This option is typically used to search a 
subspace exhaustively for vectors generating a nontrivial invariant 
subspace. 
Of course "-1" and "-m" cannot be used together. Also, "-m" cannot be 
used together with "-n".


@subsection srchmode  Specifying the Search Mode
What @b zsp does after spinning up a seed vector depends on the options
"-e", and "-c". Without any of these options, @b zsp tries to find a proper
invariant subspace. If the seed vector generates the whole space, @b zsp
tries the next seed vector and repeats until a proper invariant subspace
has been found, or until there are no more seed vectors.

With "-e", @b zsp tries to find a cyclic vector. In this mode, the program 
spins up seed vectors one-by-one until it finds a vector that generates 
the whole space, or until there are no more seed vectors available.

If you use the the option "-c" instead, @b zsp combines the span of all
seed vectors. In other words, "-c" calculates the closure of the
seed space under the generators.  For example
||    zsp -c -b sub seed gen1 gen2
calculates the closure of "seed" under the two generators and writes
a basis of the invariant subspace to "sub". 
@b zsp will print an error message if you try to use "-c"
together with any of "-1", "-m", or "-e".

Using the "-d" option you can set an upper limit on the subspace
dimension. When @b zsp finds an invariant subspace, it will stop
searching only if the dimension is at most @em MaxDim. Otherwise
the search continues with the next seed vector. "-d" cannot be used
together with neither "-e" nor "-c".


@subsection zsp_stdb Standard Basis
If you use "-t", @b zsp spins up canonically, producing the "standard 
basis". In this mode, the production of the subspace from the seed 
vector is independent of the chosen basis. Note that the standard
basis algorithm allocates an additional matrix of the same size as
the generators.


@subsection zsp_output_options Specifying Output Files
@b zsp can produce four different output files, which are all optional.
If you use the "-b" option, a basis of the invariant subspace is written
to @em Bas. The basis is always in echelon form. 

"-s" and "-q" tell @b zsp to calculate the action of the generators on 
the subspace and on the quotient, respectively. The file names are 
treated as base names with the same convention as explained above. 
For example,
<pre>
zsp -q quot -s sub gen1 gen2 seed
</pre>
Finds an invariant subspace, calculates the action on subspace an 
quotient, and write the action to "sub.1", "sub.2", "quot.1", and 
"quot.2". A second example:
<pre>
zsp -c -s std -g 3 gen pw
</pre>
Here, a standard basis is constructed using three generators, "gen.1",
"gen.2", and "gen.3", and seed vectors from "pw". The generators are 
then transformed into the standard basis and written to "std.1",
"std.2", and "std.3".

Note that "-s" and "-q" can only be used if the generators are
matrices. If you spin up with permutations, use "-b" to make a basis
of the invariant subspace, then calculate the action of the
generators using ZMU and ZCL. Example:
<pre>
zsp -b sub perm1 perm2 seed
zmu sub perm1 img
zcl sub img dummy sub1
zmu sub perm2 img
zcl sub img dummy sub2
</pre>
After these commands, sub1 and sub2 contain the action of the permutations
on the subspace.

Finally, you can write a spin-up script by using the "-o" option.
The spin-up script contains the operations performed by the spin-up
algorithm to create the subspace from the seed vectors and the 
generators. It can be used with the @ref prog_zsc "zsc" program to
repeat the same process with different seed vectors and generators.
Details on the format of the spin-up script can be found in the library reference
under SpinUp().


@section zsp_impl Implementation Details
All generators, the seed vectors (depending on "-1" and "-n"), and
a workspace are hold in memory.
The workspace is the size as generators unless the maximal
dimension has been restricted with "-d". In standard basis mode,
an additional matrix of the same size as the generators is allocated.

**/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
