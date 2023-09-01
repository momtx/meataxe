////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Condense a permutation
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <stdlib.h>
#include <string.h>


/* ------------------------------------------------------------------
   Global variables
   ------------------------------------------------------------------ */


static MtxApplicationInfo_t AppInfo = { 
"zkd", "Condense a permutation", 
"SYNTAX\n"
"    zkd [-QV] <Field> <Orbits> <Perm> <Kond>\n"
"\n"
"ARGUMENTS\n"
"    <Field> ................. The field to use for condensation\n"
"                              or 'Z' to condense over the integers.\n"
"    <Orbits> ................ Name of orbit sizes file.\n"
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
static uint32_t *RowZ;	/* One row of the ouput (for integer kondensation) */
static MtxApplication_t *App = NULL;



/* ------------------------------------------------------------------
   init() - Initialize everything
   ------------------------------------------------------------------ */


static int Init(int argc, char **argv)

{
    App = appAlloc(&AppInfo,argc,argv);
    if (App == NULL)
	return -1;

    /* Command line.
       ------------- */
    if (appGetArguments(App,4,4) < 0)
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
    if ((orbfile = sysFopen(orbname,"rb")) == NULL)
	return -1;
    Orbits = imatRead(orbfile);
    if (Orbits == NULL)
    {
	mtxAbort(MTX_HERE,"Error reading orbit table from %s",orbname);
	return -1;
    }
    OrbitSizes = imatRead(orbfile);
    if (OrbitSizes == NULL)
    {
	mtxAbort(MTX_HERE,"Error reading orbit sizes table from %s",orbname);
	return -1;
    }
    NOrbits = OrbitSizes->Noc;

    /* Read the permutation.
       --------------------- */
    Perm = permLoad(permname);
    if (Perm == NULL)
	return -1;
    if (Perm->Degree != Orbits->Noc)
    {
	mtxAbort(MTX_HERE,"%s and %s: %s\n",permname,orbname,MTX_ERR_INCOMPAT);
	    return -1;
    }
    Degree = Perm->Degree;

    /* Set field and allocate output buffer.
       ------------------------------------- */
    if (fl != 0)		/* Condensation ofer GF(q) */
    {
	ffSetField(fl);
	MESSAGE(1,("Condensation over GF(%d), characteristic is %d\n",
	    fl,ffChar));
    	m1 = ffAlloc(1, NOrbits);
    	hsz = ffAlloc(1, NOrbits);

    	/* Find the largest power of the characteristic
           which divides any of the orbit sizes
           -------------------------------------------- */
        ppow = ffChar;
    	for (i = 0; i < NOrbits; ++i)
    	{
	    n = OrbitSizes->Data[i];
	    while (n % ppow == 0)
	    	ppow *= ffChar;
    	}
    	ppow /= ffChar;
    	MESSAGE(0,("p-part taken has order %d\n",ppow));
    }
    else			/* Condensation over Z */
    {
	MESSAGE(1,("Condensation over Z\n"));
	RowZ = NALLOC(uint32_t,NOrbits);
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

    ffMulRow(hsz,FF_ZERO, NOrbits);
    for (i = 0; i < NOrbits; ++i)
    {	
	int l = (int) OrbitSizes->Data[i];
	if (l % ppow == 0)
	    f = ffInv(ffFromInt((l / ppow) % ffChar));
	else
	    f = FF_ZERO;
	ffInsert(hsz,i,f);
    }
}



int main(int argc, char **argv)
{
    int orbit;	/* Current orbit */

    if (Init(argc,argv) != 0)
    {
	mtxAbort(MTX_HERE,"Initialization failed");
	return 1;
    }
    if (readdata() != 0)
    {
	mtxAbort(MTX_HERE,"Error reading input files");
	return 1;
    }
    if (fl != 0)
    {
	init2();
    	if ((kondfile = mfCreate(kondname,fl,NOrbits,NOrbits)) == NULL)
	    return 1;
    }
    else
    {
    	if ((kondfile = mfCreate(kondname,-8,NOrbits,NOrbits)) == NULL)
	    return 1;
    }


    for (orbit = 0; orbit < NOrbits; ++orbit)
    {
	int pt;

	/* Clear the output buffer.
	   ------------------------ */
	if (fl != 0) {
	    ffMulRow(m1,FF_ZERO, NOrbits);
        }
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
		f = ffExtract(m1,image_orbit);
		f = ffAdd(f,ffExtract(hsz,image_orbit));
		ffInsert(m1,image_orbit,f);
	    }
	    else
	    {
		 RowZ[image_orbit]++;
	    }
	}

	/* Write one row to the output file.
	   --------------------------------- */
	if (fl != 0)
	    mfWriteRows(kondfile,m1,1);
	else
	    mfWrite32(kondfile, RowZ, NOrbits);
    }

    mfClose(kondfile);
    appFree(App);
    return 0;
}



/**
@page prog_zkd zkd - Condense a Permutation

@section zkd_syntax Command Line
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

@section zkd_inp Input Files
@par @em Orbits
Orbits and orbit sizes.
@par @em Perm
Permutation to be condensed.


@section zkd_out Output Files
@par @em Kond
Condensed permutation.


@section zkd_desc Description
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


@section zkd_impl Implementation Details
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
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
