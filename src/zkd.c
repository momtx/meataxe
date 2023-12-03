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
static int fieldOrder = -1;	/* Field order, 0 = integer condensation */
static int ppow;	/* l.c.m. of the p-parts of orbit sizes */
static int Degree;	/* Degree of the permutation */
static IntMatrix_t *orbits;
static IntMatrix_t *orbitSizes;
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
    if (!strcmp(App->argV[0],"Z"))
	fieldOrder = 0;
    else
	fieldOrder = atoi(App->argV[0]);
    orbname = App->argV[1];
    permname = App->argV[2];
    kondname = App->argV[3];
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int readdata()
{
   int i, n;

   // Read the orbit file.
   MtxFile_t* fileOrb = mfOpen(orbname, "rb");
   orbits = imatRead(fileOrb);
   orbitSizes = imatRead(fileOrb);
   NOrbits = orbitSizes->noc;
   mfClose(fileOrb);

   // Read the permutation.
   Perm = permLoad(permname);
   if (Perm->degree != orbits->noc)
      mtxAbort(MTX_HERE, "%s and %s: %s", permname, orbname, MTX_ERR_INCOMPAT);
   Degree = Perm->degree;

   // Set field and allocate output buffer.
   if (fieldOrder != 0) {               /* Condensation over GF(q) */
      ffSetField(fieldOrder);
      MTX_LOGD("Condensation over GF(%d), characteristic is %d", fieldOrder, ffChar);
      m1 = ffAlloc(1, NOrbits);
      hsz = ffAlloc(1, NOrbits);

      /* Find the largest power of the characteristic
         which divides any of the orbit sizes
         -------------------------------------------- */
      ppow = ffChar;
      for (i = 0; i < NOrbits; ++i) {
         n = orbitSizes->data[i];
         while (n % ppow == 0) {
            ppow *= ffChar;
         }
      }
      ppow /= ffChar;
      MTX_LOGI("p-part taken has order %d", ppow);
   } else {                     /* Condensation over Z */
      MTX_LOGD("Condensation over Z");
      RowZ = NALLOC(uint32_t, NOrbits);
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
	int l = (int) orbitSizes->data[i];
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
    if (fieldOrder != 0)
    {
	init2();
    	if ((kondfile = mfCreate(kondname,fieldOrder,NOrbits,NOrbits)) == NULL)
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
	if (fieldOrder != 0) {
	    ffMulRow(m1,FF_ZERO, NOrbits);
        }
	else
	    memset(RowZ,0,sizeof(RowZ[0]) * NOrbits);

	for (pt = 0; pt < Degree; ++pt)
	{
	    int image_orbit;
	    FEL f;
	    if (orbits->data[pt] != orbit)	    /* Ignore other orbits */
		continue;

	    /* Find out to which orbit <pt> is mapped.
	       --------------------------------------- */
	    image_orbit = orbits->data[Perm->data[pt]];

	    /* Update the condensed row.
	       ------------------------- */
	    if (fieldOrder != 0)
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
	if (fieldOrder != 0)
	    ffWriteRows(kondfile,m1,1, NOrbits);
	else
	    mfWrite32(kondfile, RowZ, NOrbits);
    }

    mfClose(kondfile);
    appFree(App);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

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
Let r be the number of orbits, O<sub>1</sub>,…,O<sub>r</sub> the orbits,
and l<sub>i</sub>:=|O_i| the size of the i-th orbit.
The first step is to calculate the largest power (m) of the
characteristic that divides any of the orbit sizes. ZKD assumes that
this is the order of the Sylow-p subgroup of the condensation
subgroup, but it prints out its findings with the message
```
p-part taken has order N
```
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

   e = 1/|K|·∑<sub>h∈K</sub>h

in FG associated to K. Now, let V be a FG-module, for example (as in this program) the natural
permutation module V=Fⁿ, where G acts by permuting the entries of vectors.
Then, Ve is an e(FG)e-module, and for any π∈G, the condensed form is eπe,
regarded as a linear map on Ve.

To calculate the action of eπe, let (v_1,…,v_n) be the standard basis such that
v<sub>i</sub>·π=v<sub>iπ</sub> for π∈G.
A basis of Ve is given by the orbit sums

   ew<sub>i</sub> = ∑<sub>k∈O<sub>i</sub></sub>v<sub>k</sub> (1≤i≤r)

and with respect to this basis we have
       
   w<sub>i</sub>(eπe) = ∑<sub>k∈O<sub>i</sub></sub> 1/(l<sub>[kπ]</sub>) w<sub>[kπ]</sub>

where [m] denotes the orbit containing m.

If K is not a p'-subgroup, e is no longer defined. However, the last formula can still
be given a sense by replacing

   1/l<sub>[kπ]</sub> → 1/l<sub>[kπ]/p<sup>m</sup></sub>  if p<sup>m</sup> | l<sub>[kπ]</sub><br>
   1/l<sub>[kπ]</sub> → 0  otherwise

where m is the highest power of the characteristic which divides any of the orbit sizes.
Thus, all but the orbits with maximal p-part are discarded, and the corresponding columns
in the output matrix are zero.
*/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
