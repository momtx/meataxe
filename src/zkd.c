/* ============================= C MeatAxe ==================================
   File:        $Id: zkd.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Condense a permutation
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <stdlib.h>
#include <string.h>


/* ------------------------------------------------------------------
   Global variables
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

static MtxApplicationInfo_t AppInfo = { 
"zkd", "Condense a permutation", 
"SYNTAX\n"
"    zkd [-QV] <Field> <Orbits> <Perm> <Kond>\n"
"\n"
"ARGUMENTS\n"
"    <Field> ................. The field to use for condensation\n"
"                              or 'Z' to condense over the integers.\n"
"    <Orbits> ................ Name of Orbit sizes file.\n"
"    <Perm> .................. Permutation to be condensed.\n"
"    <Kond> .................. File name for condensed permutation."
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"\n"
"FILES\n"
"    <Orbits> ................ I Orbit file producey by ZMO\n"
"    <Perm> .................. I Permutation to be condensed\n"
"    <Kond> .................. I Condensed permutation (square matrix)\n"
};

static const char *orbname, *permname, *kondname;
static MtxFile_t *kondfile;
static int fl = -1;	/* Field order, 0 = integer condensation */
static int ppow;	/* l.c.m. of the p-parts of orbit sizes */
static int Degree;	/* Degree of the permutation */
static IntMatrix_t *Orbits;
static IntMatrix_t *OrbitSizes;
static int NOrbits;
static PTR hsz;
static Perm_t *Perm;	/* The permutation to be condensed */
static PTR m1;		/* One row of the output matrix */
static long *RowZ;	/* One row of the ouput (for integer kondensation) */
static MtxApplication_t *App = NULL;



/* ------------------------------------------------------------------
   init() - Initialize everything
   ------------------------------------------------------------------ */


static int Init(int argc, const char **argv)

{
    App = AppAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;

    /* Command line.
       ------------- */
    if (AppGetArguments(App,4,4) < 0)
	return -1;
    if (!strcmp(App->ArgV[0],"Z"))
	fl = 0;
    else
	fl = atoi(App->ArgV[0]);
    orbname = App->ArgV[1];
    permname = App->ArgV[2];
    kondname = App->ArgV[3];
    return 0;
}




/* ------------------------------------------------------------------
   readdata() - Open and read input files
   ------------------------------------------------------------------ */

static int readdata()

{
    FILE *orbfile;
    int i, n;

    /* Read the orbit file.
       -------------------- */
    if ((orbfile = SysFopen(orbname,FM_READ)) == NULL)
	return -1;
    Orbits = ImatRead(orbfile);
    if (Orbits == NULL)
    {
	MTX_ERROR1("Error reading orbit table from %s",orbname);
	return -1;
    }
    OrbitSizes = ImatRead(orbfile);
    if (OrbitSizes == NULL)
    {
	MTX_ERROR1("Error reading orbit sizes table from %s",orbname);
	return -1;
    }
    NOrbits = OrbitSizes->Noc;

    /* Read the permutation.
       --------------------- */
    Perm = PermLoad(permname);
    if (Perm == NULL)
	return -1;
    if (Perm->Degree != Orbits->Noc)
    {
	MTX_ERROR3("%s and %s: %E\n",permname,orbname,MTX_ERR_INCOMPAT);
	    return -1;
    }
    Degree = Perm->Degree;

    /* Set field and allocate output buffer.
       ------------------------------------- */
    if (fl != 0)		/* Condensation ofer GF(q) */
    {
	FfSetField(fl);
	FfSetNoc(NOrbits);
	MESSAGE(1,("Condensation over GF(%d), characteristic is %d\n",
	    fl,FfChar));
    	m1 = FfAlloc(1);
    	hsz = FfAlloc(1);

    	/* Find the largest power of the characteristic
           which divides any of the orbit sizes
           -------------------------------------------- */
        ppow = FfChar;
    	for (i = 0; i < NOrbits; ++i)
    	{
	    n = OrbitSizes->Data[i];
	    while (n % ppow == 0)
	    	ppow *= FfChar;
    	}
    	ppow /= FfChar;
    	MESSAGE(0,("p-part taken has order %d\n",ppow));
    }
    else			/* Condensation over Z */
    {
	MESSAGE(1,("Condensation over Z\n"));
	RowZ = NALLOC(long,NOrbits);
    }

    return 0;
}


/* ------------------------------------------------------------------
   init2() - Initialize hsz, calculate starting points
   ------------------------------------------------------------------ */

static void init2()

{
    int i;
    FEL f;

    FfMulRow(hsz,FF_ZERO);
    for (i = 0; i < NOrbits; ++i)
    {	
	int l = (int) OrbitSizes->Data[i];
	if (l % ppow == 0)
	    f = FfInv(FfFromInt((l / ppow) % FfChar));
	else
	    f = FF_ZERO;
	FfInsert(hsz,i,f);
    }
}



int main(int argc, const char **argv)

{
    int orbit;	/* Current orbit */

    if (Init(argc,argv) != 0)
    {
	MTX_ERROR("Initialization failed");
	return 1;
    }
    if (readdata() != 0)
    {
	MTX_ERROR("Error reading input files");
	return 1;
    }
    if (fl != 0)
    {
	init2();
    	if ((kondfile = MfCreate(kondname,fl,NOrbits,NOrbits)) == NULL)
	    return 1;
    }
    else
    {
    	if ((kondfile = MfCreate(kondname,-8,NOrbits,NOrbits)) == NULL)
	    return 1;
    }


    for (orbit = 0; orbit < NOrbits; ++orbit)
    {
	int pt;

	/* Clear the output buffer.
	   ------------------------ */
	if (fl != 0)
	    FfMulRow(m1,FF_ZERO);
	else
	    memset(RowZ,0,sizeof(RowZ[0]) * NOrbits);

	for (pt = 0; pt < Degree; ++pt)
	{
	    int image_orbit;
	    FEL f;
	    if (Orbits->Data[pt] != orbit)	    /* Ignore other orbits */
		continue;

	    /* Find out to which orbit <pt> is mapped.
	       --------------------------------------- */
	    image_orbit = Orbits->Data[Perm->Data[pt]];

	    /* Update the condensed row.
	       ------------------------- */
	    if (fl != 0)
	    {
		f = FfExtract(m1,image_orbit);
		f = FfAdd(f,FfExtract(hsz,image_orbit));
		FfInsert(m1,image_orbit,f);
	    }
	    else
	    {
		 RowZ[image_orbit]++;
	    }
	}

	/* Write one row to the output file.
	   --------------------------------- */
	if (fl != 0)
	    MfWriteRows(kondfile,m1,1);
	else
	    MfWriteLong(kondfile,RowZ,NOrbits);
    }

    MfClose(kondfile);
    AppFree(App);
    return 0;
}



/**
@page prog_zkd zkd - Condense a Permutation

<<<<<<< HEAD
@section syntax Command Line
=======
@section zkd_syntax Command Line
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
<pre>
zkd @em Options @em Field @em Orbits @em Perm @em Kond
</pre>

@par @em Options
Standard options, see @ref prog_stdopts

@par @em Field
The field order or "Z" for integer condensation.
@par @em Orbits
Orbits and orbit sizes.
@par @em Perm
Permutation to be condensed.
@par @em Kond
Condensed permutation.

<<<<<<< HEAD
@section inp Input Files
=======
@section zkd_inp Input Files
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
@par @em Orbits
Orbits and orbit sizes.
@par @em Perm
Permutation to be condensed.


<<<<<<< HEAD
@section out Output Files
=======
@section zkd_out Output Files
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
@par @em Kond
Condensed permutation.


<<<<<<< HEAD
@section desc Description
=======
@section zkd_desc Description
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
This program reads an orbit file (@em Orbits) and a permutation
from @em Perm.
It outputs the condensed form, i.e., a matrix over GF(q) to @em Result.
The field must be specified on the command line because the other input
data is is all to do with permutations 
and the program would otherwise not know which field was intended.
The orbit file must contain two integer matrices containing the orbit 
numbers for each point and the orbit sizes, repectively. It is usually
produced by the @ref prog_zmo "zmo" program.

The second input file, @em Perm, must contain one or more permutations.
Notice that only the first permutation is read in and condensed.
If there are more than one permutation, the others are ignored.
Unlike in previous versions of this program, it is not assumed that the
orbits are contiguous.

@par Integer condensation
If @em Field is the letter "Z", @b zkd condenses over the integers.
In this case, @em Result is an integer matrix with the same dimensions as in the
GF(q) case.


<<<<<<< HEAD
@section impl Implementation Details
=======
@section zkd_impl Implementation Details
>>>>>>> 4a68ae339f0300470810ab3c90387657ccf21f0c
Let r be the number of orbits,
@f$O_1,\ldots,O_r@f$ the orbits and, for @f$i=1,\ldots,r@f$, @f$l_i:=|O_i|@f$
the size of the i-th orbit.
The first step is to calculate the largest power (m) of the
characteristic that divides any of the orbit sizes. ZKD assumes that
this is the order of the Sylow-p subgroup of the condensation
subgroup, but it prints out its findings with the message
<pre>
p-part taken has order N
</pre>
so the user can check it. If this is not the order of the Sylow-p
subgroup of the condensation group, the program will not know, so
will continue. Normally, however, the condensation subgroup K will
have trivial Sylow-p subgroup, or at any rate the Sylow subgroup
will have a regular orbit, and in this case at least the condensation
is legitimate.

The output is a square matrix with one row and one column for each
orbit of K. Abstractly, the condensation can be described as follows.
Let G be a permutation group of degree n, F a field of
characteristic p and K≤G a p'-subgroup. Then, there is an idempotent
@f[
        e = \frac{1}{|K|} \sum_{h\in K} h \in FG
@f]
associated to K. Now, let V be a FG-module, for example (as in
this program) the natural permutation module @f$V=F^n@f$, where G acts
by permuting the entries of vectors. Then, Ve is an e(FG)e-module,
and for any π∈G, the condensed form is eπe, regarded as a linear map on Ve.

To calculate the action of eπe, let @f$(v_1,\ldots,v_n)@f$ be the
standard basis such that @f$v_i\pi=v_{i\pi}@f$ for π∈G.
A basis of Ve is given by the orbit sums
@f[
        w_i = \sum_{k\in O_i} v_k
    	\qquad(1\leq i\leq r)
@f]
and with respect to this basis we have
@f[
        w_i (e\pi e) = \sum_{k\in O_i} \frac{1}{l_{[k\pi]}} w_{[k\pi]}
@f]
where [m] denotes the orbit containing m.

If K is not a p'-subgroup, e is no longer defined. However, the
last formula can still be given a sense by replacing
@f[
        \frac{1}{l_{[i\pi]}} \to \lambda_{[i\pi]}:=
        \left\{\begin{array}{ll}
        \frac{1}{l_{[i\pi]}/p^m}  & \mbox{if~}p^m|l_{[i\pi]}\\
        0                         & \mbox{otherwise}
        \end{array}\right.
@f]
where m is the highest power of the characteristic which divides
any of the orbit sizes. Thus, all but the orbits with maximal p-part
are discarded, and the corresponding columns in the output matrix are zero.
*/
