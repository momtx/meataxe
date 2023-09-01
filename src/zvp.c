////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Vector permute (make permutations from matrices)
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <string.h>
#include "meataxe.h"

#define MAXVEC 100000   /* Default value for maxvec */

#define MAX_GENERATORS  50

////////////////////////////////////////////////////////////////////////////////////////////////////
// Global data


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

static int NGen = 2;                    // number of generators
static Matrix_t *Gen[MAX_GENERATORS];   // generators
static Matrix_t *Seed;                  // seed space
static int Generate = 0;                // generate vectors (-m)
static int SeedStart = 0;               // first seed vector to use (-s)
static int maxvec = MAXVEC;             // max. orbit size (-l)

static int tabsize;                     // size of hash table
static int *vecpos;                     // position of vectors in table
static int *vecno;                      // numbers of vectors in table
static char *isfree;
static int nvec;                        // number of vectors obtained so far
static int nfinished;                   // number of finished vectors
static PTR vtable;                      // holds the vectors
static PTR tmp;
static uint32_t *perm[MAX_GENERATORS];      // permutations

static long iseed = 0;                  // current seed vector
int proj = 0;                           // operate on 1-spaces (-p)
static int vecout = 0;                  // output vectors, too (-v)
static int noout = 0;                   // no output, print orbit sizes only (-n)
static long hasm, hasp, neh;            // parameters for the hash function

static const char *MatName;             // matrix name
static const char *PermName;            // permutation name
static const char *SeedName;            // seed space name
static const char *OrbName = "orbit";   // orbit file name (optional)

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Initialize the hash parameters

static void InitHash()
{
   int dim = Gen[0]->Noc;
   int tod, i;

   neh = (dim > 25) ? 25 : dim;
   tod = 100;
   hasp = tabsize;
   if (hasp < 110) {
      tod = 7;
   }
   if (hasp > 11) {
      i = 1;
      while (++i <= tod) {
         if (hasp % i == 0) {
            --hasp;
            i = 1;
         }
      }
   }
   if (hasp < 100) {
      hasm = 3;
   } else {
      if (ffChar % 83 == 0) {
         hasm = 89;
      } else {
         hasm = 83;
      }
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static int ParseCommandLine()
{
   // options
   noout = appGetOption(App,"-n");
   vecout = appGetOption(App,"-v");
   proj = appGetOption(App,"-p");
   NGen = appGetIntOption(App,"-g",2,1,MAX_GENERATORS);
   maxvec = appGetIntOption(App,"-l",MAXVEC,0,-1);
   Generate = appGetOption(App,"-m");
   SeedStart = appGetIntOption(App,"-s",1,1,10000000) - 1;

   // arguments
   if (appGetArguments(App,3,4) < 0) {
      return -1;
   }
   MatName = App->ArgV[0];
   SeedName = App->ArgV[1];
   PermName = App->ArgV[2];
   if (App->ArgC > 3) {
      OrbName = App->ArgV[3];
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int ReadFiles()
{
   int i;

   // read generators
   for (i = 0; i < NGen; ++i) {
      char fn[200];
      sprintf(fn,"%s.%d",MatName,i + 1);
      if ((Gen[i] = matLoad(fn)) == NULL) {
         return -1;
      }
      if (Gen[i]->Nor != Gen[i]->Noc) {
         mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_NOTSQUARE);
         return -1;
      }
      if ((i > 0) && ((Gen[i]->Field != Gen[0]->Field)
                      || (Gen[i]->Nor != Gen[0]->Nor))) {
         mtxAbort(MTX_HERE,"%s and %s.%d: %s",fn,MatName,i + 1,MTX_ERR_INCOMPAT);
         return -1;
      }
   }

   // read seed space
   if ((Seed = matLoad(SeedName)) == NULL) {
      return -1;
   }
   if ((Seed->Field != Gen[0]->Field) || (Seed->Noc != Gen[0]->Nor)) {
      mtxAbort(MTX_HERE,"%s and %s.0: %s",SeedName,MatName,MTX_ERR_INCOMPAT);
      return -1;
   }

   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static int AllocateTables()
{
   int i;

   tabsize = maxvec + maxvec / 10 + 1;
   MESSAGE(1,("Allocating tables (size=%d)\n",tabsize));
   vtable = ffAlloc(tabsize + 1, Seed->Noc);
   tmp = ffAlloc(1, Seed->Noc);
   vecpos = NALLOC(int,tabsize + 1);
   vecno = NALLOC(int,tabsize + 1);
   isfree = NALLOC(char,tabsize + 1);
   for (i = 0; i < NGen; ++i) {
      perm[i] = NALLOC(uint32_t,maxvec + 1);
   }

   InitHash();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int Init(int argc, char **argv)
{
   App = appAlloc(&AppInfo,argc,argv);
   if (App == NULL) {
      return -1;
   }
   if (ParseCommandLine() != 0) {
      return -1;
   }
   if (ReadFiles() != 0) {
      return -1;
   }
   if (AllocateTables() != 0) {
      return -1;
   }
   iseed = SeedStart;
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Normalize a vector

static void Normalize(PTR row)
{
   FEL f;
   ffFindPivot(row,&f, Seed->Noc);
   ffMulRow(row,ffInv(f), Seed->Noc);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Make the next seed vector in <tmp>

static int MakeNextSeedVector()
{
   MESSAGE(1,("Starting with seed vector %ld\n",iseed));
   if (Generate) {
      iseed = MakeSeedVector(Seed,iseed,tmp);
      if (iseed < 0) {
         return -1;
      }
   } else {
      if ((int) iseed >= Seed->Nor) {
         return -1;
      }
      ffCopyRow(tmp,matGetPtr(Seed,(int)iseed), Seed->Noc);
      ++iseed;
   }
   if (proj) {
      Normalize(tmp);
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// The hash function

static int hash(PTR row)
{
   register int pos = 0, i;

   for (i = 0; i < neh; ++i) {
      pos = (pos * hasm + ffToInt(ffExtract(row,i))) % hasp;
   }
   return pos;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Prepare everything for spin-up. Assumes that the seed vector is in tmp.

static void InitTables()
{
   int i;
   int pos;
   PTR x;

   for (i = 0; i < (int) tabsize; ++i) {
      isfree[i] = 1;
   }
   nvec = 1;
   nfinished = 0;
   pos = hash(tmp);
   x = ffGetPtr(vtable,pos,Seed->Noc);
   ffCopyRow(x,tmp, Seed->Noc);
   isfree[pos] = 0;
   vecpos[0] = pos;
   vecno[pos] = 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static int MakeOrbit()
{
   int pos, im, pos1;
   PTR x;
   int igen = 0;    // generator to apply

   while (nfinished < nvec && nvec <= maxvec) {
      MESSAGE(3,("Vec[%d] * Gen[%d] = ",nfinished,igen));
      x = ffGetPtr(vtable,vecpos[nfinished],Seed->Noc);
      ffMapRow(x,Gen[igen]->Data,Seed->Noc,Seed->Noc,tmp);
      if (proj) {
         Normalize(tmp);
      }

      // Look up the vector in the hash table.
      pos1 = pos = hash(tmp);
      x = ffGetPtr(vtable,pos,Seed->Noc);
      while (!isfree[pos] && ffCmpRows(tmp,x,Seed->Noc)) {
         if (++pos == tabsize) {
            pos = 0;
            x = vtable;
         } else {
            ffStepPtr(&x,Seed->Noc);
         }
         // always true since the hash table is always larger than maxvec
         MTX_ASSERT(pos != pos1);
      }

      if (isfree[pos]) {        // new vector
         MESSAGE(3,("%d (new)\n",nvec));
         ffCopyRow(x,tmp, Seed->Noc);
         im = nvec;
         isfree[pos] = 0;
         vecpos[nvec] = pos;
         vecno[pos] = nvec;
         nvec++;
      } else {
         im = vecno[pos];
         MESSAGE(3,("%d\n",im));
      }
      if (nvec % 10000 == 0) {
         MESSAGE(2,("%d vectors, %d finished\n",nvec,nfinished));
      }
      perm[igen][nfinished] = im;

      // next generator
      if (++igen >= NGen) {
         igen = 0;
         ++nfinished;
      }
   }

   if (nfinished < nvec) {
      return -1;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void Cleanup()
{
   appFree(App);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static int WriteOutput()
{
   int rc = 0;
   int i;
   char fn[200];

   if (noout) {
      return 0;
   }

   /* Write vectors
      ------------- */
   if (vecout) {
      MtxFile_t *f;
      int i;

      MESSAGE(1,("Writing orbit to %s\n",OrbName));
      if ((f = mfCreate(OrbName,ffOrder,nvec,Seed->Noc)) == 0) {
         mtxAbort(MTX_HERE,"Cannot open %s",OrbName);
         rc = -1;
      }
      for (i = 0; i < nvec; ++i) {
         PTR x = ffGetPtr(vtable,vecpos[i],Seed->Noc);
         mfWriteRows(f,x,1);
      }
      mfClose(f);
   }

   /* Write permutations
      ------------------ */
   MESSAGE(1,("Writing permutations\n"));
   for (i = 0; i < NGen; ++i) {
      sprintf(fn,"%s.%d",PermName,i + 1);
      MtxFile_t* f = mfCreate(fn,MTX_TYPE_PERMUTATION, nvec, 1);
      mfWrite32(f,perm[i],nvec);
      mfClose(f);
   }
   return rc;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   int rc = 1;

   if (Init(argc,argv) != 0) {
      mtxAbort(MTX_HERE,"Initialization failed");
      return 1;
   }

   /* Main loop
      --------- */
   while (1) {
      if (MakeNextSeedVector() != 0) {
         mtxAbort(MTX_HERE,"No nore seed vectors");
         break;
      }
      InitTables();
      if (MakeOrbit() == 0) {
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

@section zvp_syntax Command Line
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

@section zvp_inp Input Files

@par @em Mat.1, @em Mat.2, ...
Generators (square matrices).

@par @em Seed
Seed vectors (matrix).

@section zvp_out Output Files

@par @em Perm.1, @em Perm.2, ...
Permutations.

@par @em Orbit
The orbit (matrix).

@section zvp_desc Description
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

@subsection zvp_impl Implementation Details
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
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
