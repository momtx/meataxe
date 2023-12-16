////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Calculate the dotted-lines.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <string.h>
#include <stdlib.h>

MatRep_t *Rep;			// Generators of the current constituent
Matrix_t *cycl = NULL;		// List of cyclic submodules
long *class[MAXCYCL];		// Classes of vectors
uint32_t nmountains = 0;	// Number of mountains
Matrix_t *mountlist[MAXCYCL];	// Mountains
BitString_t *subof[MAXCYCL];	// Incidence matrix
int cfstart[LAT_MAXCF+1];	// First mountain of each c.f.
char lck[MAXCYCL];
char lck2[MAXCYCL];
BitString_t *dotl[MAXDOTL];	    // Dotted lines
BitString_t *MaxMountains[MAXDOTL]; // Maximal mountains in dotted lines
int ndotl = 0;			// Number of dotted-lines in <dotl>
int firstdotl = 0;		// Used for locking
int firstm, nextm;		// First and last+1 mountain for the
				// current constituent

// sumdim[i][j] contains the dimension of mountain[i] + mountain[j].
// Most of the program's work constists in comparing sums of pairs of
// mountains and looking if they are equal. Thus, a given sum may be 
// needed many times. Since we don't keep the sums in memory we have
// to recalculate them each time. But comparing the dimensions we can 
// avoid many unnecessary calculations.
long *sumdim[MAXCYCL];

// This is the length of a dotted line for the current constituent.
// For theoretical reasions the is always Q+1, where GF(Q) is the
// splitting field for the constituent.
int dotlen;

static int writeGapOutput = 0;
static int Opt_FindDuplicates = 0;  // Find 'duplicate' dotted-lines
static Lat_Info LI;		    // Data from .cfinfo file



static MtxApplicationInfo_t AppInfo = { 
"mkdotl", "Find Dotted-Lines", 
"\n"
"SYNTAX\n"
"    mkdotl [<Options>] <Name>\n"
"\n"
"ARGUMENTS\n"
"    <Name> .................. Name of the representation\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -G ...................... GAP output (implies -Q)\n"
"    --nodup ................. Find and discard duplicate dotted lines\n"
"\n"
"FILES\n"
"    <Name>.cfinfo ........... IO Constituent info file\n"
"    <Name><Cf>.v ............ I  Cyclic submodules, generated by MKCYCL\n"
"    <Name>.inc .............. I  Incidence matrix generated by MKINC\n"
"    <Name>.mnt .............. I  Mountain data (from MKINC)\n"
"    <Name>.dot .............. O  Dotted-lines\n"
};

static MtxApplication_t *App = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////

static void readFiles(const char* basename)
{
   char fn[200];
   int i;

   latReadInfo(&LI, basename);

   cfstart[0] = 0;
   for (i = 1; i <= LI.nCf; ++i) {
      cfstart[i] = cfstart[i - 1] + (int)(LI.Cf[i - 1].nmount);
   }

   // Read the incidence matrix
   sprintf(fn, "%s.inc", LI.BaseName);
   MtxFile_t* f = mfOpen(fn, "rb");
   mfRead32(f, &nmountains, 1);
   if (nmountains != cfstart[LI.nCf]) {
      mtxAbort(MTX_HERE, "Bad number of mountains in .inc file");
   }
   MTX_LOGD("Reading incidence matrix (%lu mountains)", (unsigned long) nmountains);
   fflush(stdout);
   for (i = 0; i < (int) nmountains; ++i) {
      subof[i] = bsRead(f);
   }
   mfClose(f);

   for (i = 0; i < (int) nmountains; ++i) {
      sumdim[i] = NALLOC(long, nmountains);
      memset(sumdim[i], 0, (size_t)nmountains * sizeof(long));
   }

   // Read classes
   {
      sprintf(fn, "%s.mnt", LI.BaseName);
      MTX_LOGD("Reading classes (%s)", fn);
      FILE *f = sysFopen(fn, "r");
      for (i = 0; i < nmountains; ++i) {
         long mno, mdim, nvec, * p;
         int k;
         if (fscanf(f, "%ld%ld%ld", &mno, &mdim, &nvec) != 3 || mno != i || nvec < 1 || mdim < 1) {
            mtxAbort(MTX_HERE, "Invalid .mnt file");
         }
         p = class[i] = NALLOC(long, nvec + 2);
         *p++ = nvec;
         for (k = 0; k < nvec; ++k, ++p) {
            if (fscanf(f, "%ld", p) != 1 || *p < 1) {
               mtxAbort(MTX_HERE, "Invalid .mnt file");
            }
         }
         if (fscanf(f, "%ld", p) != 1 || *p != -1) {
            mtxAbort(MTX_HERE, "Invalid .mnt file");
         }
      }
      fclose(f);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void mkmount(int i)
{
    Matrix_t *seed;
    PTR x, y;
    long *p;

    seed = matAlloc(cycl->field,class[i][0],cycl->noc);
    x = seed->data;
    for (p = class[i] + 1; *p > 0; ++p)
    {
	if (*p < 1 || *p > cycl->nor)
	{
	    mtxAbort(MTX_HERE,"BAD VECTOR IN CLASS");
	    return;
	}
	y = matGetPtr(cycl,*p - 1);
	ffCopyRow(x,y, cycl->noc);
	ffStepPtr(&x, cycl->noc);
    }

    mountlist[i] = SpinUp(seed,Rep,SF_EACH|SF_COMBINE,NULL,NULL);
    if (mountlist[i] == NULL)
    {
	mtxAbort(MTX_HERE,"Cannot spin up mountain");
	return;
    }
    matFree(seed);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize everything for the next composition factor:
/// Read generators and vectors, calculate mountains...

static void initCf(int cf)
{
    char fn[200];
    int j;

    // Read the generators of the condensed module
    sprintf(fn,"%s%s.%%dk",LI.BaseName,latCfName(&LI,cf));
    Rep = mrLoad(fn,LI.NGen);

    // Read generating vectors for the cyclic submodules
    sprintf(fn,"%s%s.v",LI.BaseName,latCfName(&LI,cf));
    cycl = matLoad(fn);

    // Calculate the length of dotted-lines. This is always 
    // Q + 1 where Q is the splitting field order.
    dotlen = ffOrder;
    for (j = LI.Cf[cf].spl - 1; j > 0; --j)
    	dotlen *= ffOrder;
    ++dotlen;
    MTX_LOGD("Length of dotted-lines is %d",dotlen);

    // Calculate the mountains
    for (j = cfstart[cf]; j < cfstart[cf+1]; ++j)
	mkmount(j);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void cleanupCf()
{
    matFree(cycl);
    mrFree(Rep);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static Matrix_t *sum(int i, int k)

{
    Matrix_t *s;
    int dim_i, dim_k;

    dim_i = mountlist[i]->nor;
    dim_k = mountlist[k]->nor;

    s = matAlloc(ffOrder,dim_i + dim_k,mountlist[i]->noc);

    matCopyRegion(s,0,0,mountlist[i],0,0,dim_i,-1);
    matCopyRegion(s,dim_i,0,mountlist[k],0,0,dim_k,-1);

    matEchelonize(s);

    // Remember the dimension of the sum. We use this information 
    // later to avoid unnecessary recalculations.
    sumdim[i][k] = sumdim[k][i] = s->nor;
    return s;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void lock(int i, char *c)
{	
    int l, m;
    BitString_t *b;

    memset(c,0,sizeof(lck));
    for (m = firstm; m < nextm; ++m)
    {	
	if (bsTest(subof[i],m) || bsTest(subof[m],i))
	    c[m] = 1;
    }
    for (l = firstdotl; l < ndotl; ++l)
    {	
	b = dotl[l];
	if (!bsTest(b,i))
	    continue;
	for (m = firstm; m < nextm; ++m)
	{
	    if (bsTest(b,m))
		c[m] = 1;
	}
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void FindMaxMountains(Matrix_t *span, BitString_t *bs)
{
    int m;

    bsClearAll(bs);
    for (m = firstm; m < nextm; ++m)
    {
	Matrix_t *tmp = matDup(mountlist[m]);
	matClean(tmp,span);
	if (tmp->nor == 0) // Mountain is countained in <span>
	    bsSet(bs,m);
	matFree(tmp);
    }

    // Remove non-maximal mountains
    for (m = firstm; m < nextm; ++m)
    {
	int k;
	if (!bsTest(bs,m))
	    continue;
	for (k = firstm; k < nextm; ++k)
	{
	    if (k != m && bsTest(subof[k],m))
		bsClear(bs,k);
	}
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

// Find out if mountains #i and #k generate a dotted line.

static void trydot(int i, int k, int beg, int next)
{
    Matrix_t *span, *sp;
    int count, l, m;

/*
     OPTIMIERUNG, FUNKTIONIERT NICHT!
    if (i == 1 && k == 13)
	i = 1;

    for (l = 0; l < ndotl; ++l)
    {
	if (bsTest(MaxMountains[l],i) && bsTest(MaxMountains[l],k))
	    return;
    }
*/
    lock(k,lck2);
    BitString_t *dot = bsAlloc(nmountains);
    bsSet(dot,i);
    bsSet(dot,k);
    span = sum(i,k);
    count = 2;
    for (l = beg; l < next && count < dotlen; ++l)
    {
	int abort = 0;
	if (lck[l] || lck2[l]) 
	    continue;
	for (m = i; !abort && m < l; ++m)
	{
	    if (!bsTest(dot,m)) 
		continue;
	    if (sumdim[l][m] != 0 && sumdim[l][m] != span->nor)
	    	abort = 1;
	    else
	    {	
		sp = sum(l,m);
		abort = (sp->nor != span->nor) || !IsSubspace(span,sp,0);
		matFree(sp);
	    }
	}
	if (!abort)
	{	
	    bsSet(dot,l);
	    ++count;
	    lck[l] = 1;
	}
    }
    if (count == dotlen)	/* We have found a dotted line */
    {	
	int d;

	MTX_LOGD("New dotted line: %d+%d",i,k);
	if (ndotl >= MAXDOTL)
	{
	    mtxAbort(MTX_HERE,"Too many dotted lines (max %d)",MAXDOTL);
	    return;
	}
	dotl[ndotl] = dot;
	if (Opt_FindDuplicates)
	{
	    MaxMountains[ndotl] = bsDup(dot);
	    FindMaxMountains(span,MaxMountains[ndotl]);
	    for (d = 0; d < ndotl; ++d)
	    {
		if (bsCompare(MaxMountains[ndotl],MaxMountains[d]) == 0)
		    break;
	    }
	    if (d < ndotl)
	    {
		bsFree(dot);
		bsFree(MaxMountains[ndotl]);
		MTX_LOG2("Discarding %d+%d (= dl %d)",i,k,d);
	    }
	    else
		ndotl++;
	}
	else
	    ndotl++;
	
    }
    else
	bsFree(dot);
    matFree(span);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Find dotted-lines in one constituent

static void mkdot(int cf)
{
    int i, k;

    firstm = cfstart[cf];
    nextm = cfstart[cf+1];
    firstdotl = ndotl;
    for (i = firstm; i < nextm; ++i)
    {
	MTX_LOG2("Trying mountain %d",i);
	lock(i,lck);
    	for (k = i+1; k < nextm; ++k)
    	{	
	    if (lck[k]) 
		continue;
    	    trydot(i,k,k+1,nextm);
    	}
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Write dotted lines.

static void WriteResult()
{
   char fn[200];

   strcat(strcpy(fn, LI.BaseName), ".dot");
   MTX_LOGD("Writing %s (%d dotted line%s)", fn, ndotl, ndotl != 1 ? "s" : "");
   MtxFile_t* f = mfOpen(fn, "wb");
   const uint32_t l = ndotl;
   mfWrite32(f, &l, 1);
   for (int i = 0; i < ndotl; ++i) {
      bsWrite(dotl[i], f);
   }
   mfClose(f);
   latWriteInfo(&LI);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Write dotted lines in GAP format.

static void WriteResultGAP()
{
    int i,j ;

    printf ("MeatAxe.DottedLines := [\n") ;
    for (i = 0; i < ndotl; ++i)
    {
	printf( "BlistList([" );
	for (j = 0 ; j < nmountains ; j++)
        printf( j < (nmountains - 1) ? "%s," : "%s], [1])" , 
	       bsTest(dotl[i],j) ? "1" : "0" );
        printf( ",\n" );
    }
    printf( "];\n" );
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char** argv)
{
   App = appAlloc(&AppInfo, argc, argv);
   writeGapOutput = appGetOption(App, "-G --gap");
   Opt_FindDuplicates = appGetOption(App, "--nodup");
   appGetArguments(App, 1, 1);
   MTX_LOGI("Start mkdotl - Find dotted-lines");

   readFiles(App->argV[0]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void cleanup()
{
   latCleanup(&LI);
   for (int i = 0; i < nmountains; ++i) {
      matFree(mountlist[i]);
      bsFree(subof[i]);
   }
   for (int i = 0; i < ndotl; ++i) {
      bsFree(dotl[i]);
      if (Opt_FindDuplicates) bsFree(MaxMountains[i]);
   }
   appFree(App);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    int i, nn = 0;

    init(argc,argv);
    for (i = 0; i < LI.nCf; ++i)
    {
	initCf(i);
	mkdot(i);
	LI.Cf[i].ndotl = ndotl - nn;
	MTX_LOGI("%s%s: %d vectors, %ld mountains, %ld dotted line%s",
	    LI.BaseName,latCfName(&LI,i),  cycl->nor,LI.Cf[i].nmount,
	    LI.Cf[i].ndotl, LI.Cf[i].ndotl != 1 ? "s": "");
	nn = ndotl;
	cleanupCf();
    }

/*
    for (i = 0; i < ndotl; ++i)
    {
	int k;
	printf("%3d:",i);
	for (k = 0; k < nmountains; ++k)
	    putc(bsTest(dotl[i],k) ? 'x' : '.',stdout);
	printf("\n");
	printf("    ");
	for (k = 0; k < nmountains; ++k)
	    putc(bsTest(MaxMountains[i],k) ? 'x' : '.',stdout);
	printf("\n");
    }
*/
    WriteResult();
    if (writeGapOutput)
         WriteResultGAP();
    cleanup();
    return 0;
}


/**
@page prog_mkdotl mkdotl - Find Dotted-lines

@section mkdotl_syntax Command Line
<pre>
mkdotl [@em Options] [-G] [--nodup] @em Name
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -G
  Produce output in GAP format.
@par --nodup
  Eliminate redundant dotted-lines.
@par @em Name
  Name of the representation.

@section mkdotl_inp Input Files
@par @em Name.cfinfo
  Constituent info file.
@par @em NameCf.v
  Cyclic submodules, generated by @ref prog_mkcycl "mkcycl".
@par @em Name.inc
  Incidence matrix, created by @ref prog_mkinc "mkinc".
@par @em Name.mnt
  Mountain data, created by @ref prog_mkinc "mkinc".
  

@section mkdotl_out Output Files
@par @em Name.dot
  Dotted-lines.

@section mkdotl_desc Description
The MKDOTL program is part of the @ref sec_progs_lattice "Submodule Lattice Package".
MKDOTL calculates a set of "dotted-lines" (see @ref LMR94) between the local submodules.
More precisely, it computes one dotted-line for each submodule with head isomorphic to S⊕S,
where S is irreducible.
It can be shown that this set of dotted-lines is sufficient to determine the complete
submodule lattice as described by Benson and Conway (@ref BC85).

Input for this program are the incidence matrix calculated by
@ref prog_mkinc "mkinc" and the cyclic submodules from @ref prog_mkcycl "mkcycl".
Again, the whole calculation takes place in the condensed
modules, so there is no need to uncondense the cyclic submodules.

It is known that all dotted lines have length q+1, where q is the
order of the splitting field. This information is used by the program
to determine if a dotted line is complete.

A list of all dotted lines is written to @em Name.dot.

Using the option "--nodup" eliminates redundant dotted-lines from the output.
If this option is specified, the program will calculate, for each dotted-line,
the maximal mountains contained in the span of the dotted-line. 
If a dotted-line has the same set of maximal mountains as an earlier 
dotted-line, it is considered as redundant and dropped. Note that 
"--nodup" increases both memory and CPU time usage. However, the 
subsequent step, @ref prog_mkgraph "mkgraph",
will benefit from a reduction of the number of dotted-lines.

**/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
