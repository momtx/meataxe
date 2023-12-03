////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Find peak words and condense
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>
#include <stdlib.h>

#define MAX_MODULES 50
#define MAXCF (3 * LAT_MAXCF)       // Max. number of different constituents

int NumMods = 0;                    // Number of modules

/// @private
static struct module_struct {
   Lat_Info Info;                   // Data from .cfinfo file
   MatRep_t *Rep;                   // Generators
   WgData_t *Wg;                    // Word generators
   Matrix_t *SsBasis;               // Semisimplicity basis
} ModList[MAX_MODULES];             // List of all modules


/// Table of constituents
/// @private
struct cf_struct {
   char displayName[30];            // For messages
   MatRep_t *Gen;                   // Generators
   CfInfo *Info;                    // Constituent info
   WgData_t *Wg;                    // Word generators
   Lat_Info *Mod;                   // Which module does this come from
   int CfNo;                        // Constituent index within that module
   Matrix_t *PWNullSpace;           // Peak word null space
   int mult;                        // In how many modules does it appear?
   int CfMap[MAX_MODULES][2];       // (module,cf)
};
int NumCf = 0;                      // Number of inequivalent constituents
static struct cf_struct CfList[MAXCF];  // List of inequivalent constituents

static int opt_G = 0;
static int opt_n = 0;                 // No condensation, peak words only
static int opt_p = 0;                 // Use full polynomials in PW search
static int opt_t = 0;                 // Transform generators to std basis
static int opt_b = 0;                 // Calculate a semisimplicity basis
static int opt_k = 0;                 // Make kernel of PW

#define MAXLOCK 100
static long include[MAXLOCK][2];
static int ninclude = 0;
static long exclude[MAXLOCK][2];
static int nexclude = 0;
static int PeakWordsMissing;            // Number of missing peak words

/// @private
struct SharedData {
   // Constant data, not protected

   // Variable data, protected by mutex
   #if defined(MTX_DEFAULT_THREADS)
   pthread_mutex_t mutex;
   #endif
};
#if defined(MTX_DEFAULT_THREADS)
static struct SharedData sharedData;
#endif

static MtxApplicationInfo_t AppInfo = {
   "pwkond", "Peakword Condensation",
   "\n"
   "SYNTAX\n"
   "    pwkond [<Options>] <Name> [<Name> ...]\n"
   "\n"
   "ARGUMENTS\n"
   "    <Name> .................. Name of the representation\n"
   "\n"
   "OPTIONS\n"
   MTX_COMMON_OPTIONS_DESCRIPTION
   "    -G ...................... GAP output (implies -Q)\n"
   "    -n ...................... Find peak words only, do not condense\n"
   "    -p ...................... Use full polynomials in peak word search\n"
   "    -i <List> ............... Words to try first. Example: -i 100,20-35.\n"
   "    -e <List> ............... Exclude words from search. Example: -e 3,20-99.\n"
   "    -t ...................... Transform generators into standard basis\n"
   "    -b ...................... Calculate a semisimplicity basis\n"
   "    -k ...................... Compute kernel of peak words\n"

   "\n"
   "FILES\n"
   "    <Name>.cfinfo ........... IO Constituent info file\n"
   "    <Name>.{1,2,...} ........ I  Generators\n"
   "    <Name><Cf>.{1,2...} ..... I  Generators on the constituents\n"
   "    <Name><Cf>.{1,2...}k .... O  Condensed generators\n"
   "    <Name><Cf>.{1,2...}.std   O  Condensed generators in std basis (with -t)\n"
   "    <Name><Cf>.op ........... O  Spin-up script for standard basis (with -t)\n"
   "    <Name><Cf>.np ........... O  Condensed peak word\n"
   "    <Name><Cf>.im ........... O  Image used for condensation\n"
   "    <Name><Cf>.k ............ O  Peakword kernel (with -k or without -n)\n"
   "    <Name>.ssb .............. O  Semisimplicity basis (with -b)\n"
};

static MtxApplication_t *App = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////

// Add a constituent
// This function checks if a given constituent is already in <CfList>. If
// not, it is added to the list. If it is already in the list, <cf> is
// deleted.

static void AddConstituent(MatRep_t *cf, CfInfo *info, int modno, int cfno)
{
   int i, m;
   for (i = 0; i < NumCf; ++i) {
      if (IsIsomorphic(CfList[i].Gen,CfList[i].Info,cf,NULL,0))
         break;
   }

   if (i < NumCf) { // Constituent was already in the list
      int m = CfList[i].mult;
      mrFree(cf);
      CfList[i].CfMap[m][0] = modno;
   } else {         // It's a new constituent
      CfList[i].Gen = cf;
      CfList[i].Info = info;
      CfList[i].Wg = wgAlloc(cf);
      CfList[i].Mod = &ModList[modno].Info;
      CfList[i].mult = 0;

      ++NumCf;
   }
   m = CfList[i].mult;
   CfList[i].CfMap[m][0] = modno;
   CfList[i].CfMap[m][1] = cfno;
   CfList[i].mult++;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Adds all constituents of a module.
/// This function calls <AddConstituent()> for each constituent of the <mod>-th module in
/// <ModList>. I.e., it adds all constituents of that module to the global constituent list and
/// sets up the constituent map.

static void AddConstituents(int mod)
{
   Lat_Info *li = &ModList[mod].Info;
   int i;
   for (i = 0; i < li->nCf; ++i) {
      MatRep_t *cf;
      char* const fn = strMprintf("%s%s",li->BaseName,latCfName(li,i));
      cf = mrLoad(fn,li->NGen);
      sysFree(fn);
      AddConstituent(cf,li->Cf + i,mod,i);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void loadConstituents()
{
   for (int i = 0; i < NumMods; ++i) {
      AddConstituents(i);
   }

   // Sort the constituents by dimension to speed up the peak word search.
   for (int i = 0; i < NumCf; ++i) {
      for (int k = i + 1; k < NumCf; ++k) {
         if (CfList[i].Info->dim > CfList[k].Info->dim) {
            struct cf_struct tmp = CfList[i];
            CfList[i] = CfList[k];
            CfList[k] = tmp;
         }
      }
      struct cf_struct* const cfi = CfList + i;
      snprintf(cfi->displayName, sizeof(cfi->displayName), "cf%lu", (unsigned long)i);
   }

   for (int i = 0; i < NumCf; ++i) {
      struct cf_struct *const cf = CfList + i;
      MTX_XLOGD(sb) {
         sbPrintf(sb, "%s is", cf->displayName);
         for (size_t k = 0; k < cf->mult; ++k) {
            struct module_struct *mod = ModList + cf->CfMap[k][0];
            const Lat_Info* const li = &mod->Info;
            sbPrintf(sb, " %s%s", li->BaseName,latCfName(li,cf->CfMap[k][1]));
         }
      }
   }

}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void checkCompatibility(int i)
{
   const Lat_Info* infoI = &ModList[i].Info;
   const Lat_Info* info0 = &ModList[0].Info;
   if (infoI->NGen != info0->NGen || infoI->field != info0->field) {
      mtxAbort(MTX_HERE,"%s and %s: %s",App->argV[0],App->argV[i], MTX_ERR_INCOMPAT);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Loads the .cfinfo files and generators of all modules.

static void loadModules()
{
   // Set the number of modules.
   NumMods = App->argC;
   if (NumMods > MAX_MODULES) {
      mtxAbort(MTX_HERE,"Too many modules (max. %d allowed)",MAX_MODULES);
   }

   // Read the .cfinfo files and load the generators (if needed).
   for (int i = 0; i < NumMods; ++i) {
      latReadInfo(&ModList[i].Info,App->argV[i]);
      MTX_LOGI("%s: %d composition factors",App->argV[i],ModList[i].Info.nCf);
      checkCompatibility(i);

      // Clear any existing peak words.
      for (int k = 0; k < ModList[i].Info.nCf; ++k) {
         ModList[i].Info.Cf[k].peakWord = -1;
      }

      // Read the generators, set up ss bases and word generators.
      if (!opt_n || opt_k || opt_b) {
         ModList[i].Rep = mrLoad(App->argV[i],ModList[i].Info.NGen);
         ModList[i].Wg = wgAlloc(ModList[i].Rep);
         if (opt_b) {
            int dim = ModList[i].Rep->Gen[0]->nor;
            ModList[i].SsBasis = matAlloc(ffOrder,dim,dim);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Generalized condensation of one matrix.

static void gkond(const Lat_Info *li, int i, Matrix_t *b, Matrix_t *k, Matrix_t *w, const char *name)
{
   Matrix_t* x1 = matDup(k);
   matMul(x1,w);
   Matrix_t *x2 = QProjection(b,x1);
   char *fn = strMprintf("%s%s.%s",li->BaseName,latCfName(li,i),name);
   matSave(x2, fn);
   sysFree(fn);
   matFree(x2);
   matFree(x1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Transforms a constituent to standard basis. Writes the generators to XXX.std.N and the standard
/// basis spinup script to XXX.op.
/// Note: the generators in CfList[cf] remain unchanged!

static void transformToStandardBasis(struct cf_struct* cf)
{
   // Make the standard basis and spinup script. Transform the generators.
   MTX_LOGD("%s Transforming to standard basis", cf->displayName);
   IntMatrix_t *script = NULL;
   Matrix_t* sb =
      SpinUp(cf->PWNullSpace,cf->Gen, SF_FIRST | SF_CYCLIC | SF_STD,&script,NULL);
   MatRep_t* stdRep = mrChangeBasis2(cf->Gen, sb);
   matFree(sb);

   // Write the transformed generators and the spin-up script.
   for (int m = 0; m < cf->mult; ++m) {
      char fn[200];
      Lat_Info *li = &ModList[cf->CfMap[m][0]].Info;
      int i = cf->CfMap[m][1];
      sprintf(fn,"%s%s.op",li->BaseName,latCfName(li,i));
      MTX_LOG2("%s wrote operations to %s", cf->displayName, fn);
      imatSave(script,fn);
      for (int k = 0; k < li->NGen; ++k) {
         sprintf(fn,"%s%s.std.%d",li->BaseName,latCfName(li,i),k + 1);
         matSave(stdRep->Gen[k],fn);
      }
      MTX_LOG2("%s wrote %s%s.op and %s%s.std.(1..%d)",
            cf->displayName,
            li->BaseName,latCfName(li,i),
            li->BaseName,latCfName(li,i), li->NGen);
   }

   mrFree(stdRep);
   imatFree(script);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Find the starting point for a constituent.

static int CfPosition(const Lat_Info *li, int cf)
{
   int pos = 0;
   int i;
   MTX_ASSERT(cf >= 0 && cf < li->nCf);
   for (i = 0; i < cf; ++i) {
      pos += li->Cf[i].dim * li->Cf[i].mult;
   }
   return pos;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

// Generalized condensation for one irreducible

static void kond(struct cf_struct* cfData, struct module_struct* mod, int cf)
{
   const Lat_Info * const li = &mod->Info;
   char fn[LAT_MAXBASENAME + 10];
   Matrix_t *peakword, *kern, *m, *k, *pw;
   int j, pwr;

   // Make the peak word, find its stable power, and calculate both kernel and image.
   peakword = wgMakeWord2(mod->Wg,li->Cf[cf].peakWord);
   matInsert_(peakword,li->Cf[cf].peakPol);
   pw = matDup(peakword);
   StablePower_(peakword,&pwr,&kern);
   MTX_LOGD("%s stablePwr=%d, nul=%lu, mult=%lu, spl=%lu",
            cfData->displayName,
            pwr, (unsigned long) kern->nor,
            (unsigned long) li->Cf[cf].mult, (unsigned long) li->Cf[cf].spl);

   if (kern->nor != li->Cf[cf].mult * li->Cf[cf].spl) {
      mtxAbort(MTX_HERE,"%s Something is wrong here", cfData->displayName);
   }
   matEchelonize(peakword);

   // Write out the image
   if (!opt_n) {
      sprintf(fn,"%s%s.im",li->BaseName,latCfName(li,cf));
      matSave(peakword,fn);
   }

   // Write out the `uncondense matrix'
   m = QProjection(peakword,kern);
   k = matInverse(m);
   matFree(m);
   matMul(k,kern);
   sprintf(fn,"%s%s.k",li->BaseName,latCfName(li,cf));
   matSave(k,fn);

   // Condense all generators
   for (j = 0; j < li->NGen; ++j) {
      sprintf(fn,"%dk",j + 1);
      gkond(li,cf,peakword,k,mod->Rep->Gen[j],fn);
   }

   // Condense the peak word
   gkond(li,cf,peakword,k,pw,"np");

   // Calculate the semisimplicity basis.
   if (opt_b) {
      Matrix_t *seed, *partbas;
      int pos = CfPosition(li,cf);
      seed = matNullSpace_(pw,0);
      partbas = SpinUp(seed,mod->Rep,SF_EACH | SF_COMBINE | SF_STD,NULL,NULL);
      matFree(seed);

      if (pos < 0 || pos + partbas->nor > mod->SsBasis->nor) {
         mtxAbort(MTX_HERE,"Error making basis - '%s' is probably not semisimple", li->BaseName);
      }
      matCopyRegion(mod->SsBasis,pos,0,partbas,0,0,-1,-1);
      matFree(partbas);
   }
   matFree(pw);
   matFree(k);
   matFree(kern);
   matFree(peakword);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void condense(struct cf_struct* cf)
{
   for (int k = 0; k < cf->mult; ++k) {
      struct module_struct* const mod = ModList + cf->CfMap[k][0];
      const int i = cf->CfMap[k][1];
      MTX_LOGD("%s condensing %s%s", cf->displayName, mod->Info.BaseName, latCfName(&mod->Info,i));
      kond(cf, mod, i);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes the .cfinfo file and the semisimplicity basis.
/// It is called after a new peak word has been found.

static void WriteOutput(int final)
{
   int i;
   for (i = 0; i < NumMods; ++i) {
      latWriteInfo(&ModList[i].Info);
      MTX_LOGD("Wrote %s.cfinfo", ModList[i].Info.BaseName);
      if (opt_b) {
         char* const fn = strMprintf("%s.ssb",ModList[i].Info.BaseName);
         matSave(ModList[i].SsBasis,fn);
         MTX_LOGD("Wrote %s", fn);
         sysFree(fn);
      }
   }
   if (!final) {
      return;
   }

   // Write GAP output.
   if (opt_G) {
      int m, i;
      printf("MeatAxe.PeakWords := [\n");
      for (m = 0; m < NumMods; ++m) {
         const Lat_Info * const ModInfo = &ModList[m].Info;
         printf("# module: %s\n[\n", ModInfo->BaseName);
         for (i = 0; i < ModInfo->nCf; ++i) {
            const CfInfo * const Cf = ModInfo->Cf + i;
            printf("    # irreducible factor: %s\n", latCfName(ModInfo,i));

            StrBuffer* sb = sbAlloc(100);
            sbPrintf(sb, "    [ %ld, ", (long) Cf->peakWord);
            gapFormatWord(sb, CfList[i].Wg,Cf->peakWord);
            sbAppend(sb, ", ");
            gapFormatPoly(sb, Cf->peakPol);
            sbPrintf(sb, " ]%s\n",i == ModInfo->nCf - 1 ? "" : ",");
            sbFree(sb);
            sysFree(sb);
         }
         printf(m < NumMods - 1 ? "],\n" : "]\n");
      }
      printf("];\n");
   }

}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Copies the peak word and polynomial just found to all modules having an appropriate constituent.
/// Alo prints a nice message showing the peak word and to which constituents it applies.

static void CopyPeakWordToAllModules(struct cf_struct* cf)
{
   CfInfo * const info = cf->Info;
   const unsigned long pw = info->peakWord;
   const Poly_t* const pp = info->peakPol;

   // Copy peak word and peak polynomial to the other modules
   for (size_t k = 1; k < cf->mult; ++k) {
      CfInfo* const other = ModList[cf->CfMap[k][0]].Info.Cf + cf->CfMap[k][1];
      other->peakWord = pw;
      other->peakPol = polDup(pp);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Handle new peak word.
///
/// This function is called each time a peak word is found. Depending on the command line options
/// we condense the peak word, and transform the generators to standard basis.

// TODO:
// We also write the .cfinfo file each time. This allows the user to kill the running program and
// continue with the peak words found so far.

static void peakWordFound(void* arg)
{
   struct cf_struct*cf = (struct cf_struct*)arg;
   CfInfo* const info = cf->Info;
   const unsigned long pw = info->peakWord;

   MTX_XLOGI(sb) {
      sbPrintf(sb, "%s peakWord=%ld(%s)", cf->displayName, pw, wgSymbolicName(cf->Wg,pw));
      sbPrintf(sb, " peakPol=");
      polFormat(sb,info->peakPol);
   }

   CopyPeakWordToAllModules(cf);
   if (!opt_n || opt_k) {
      condense(cf);
   }
   if (opt_t) {
      transformToStandardBasis(cf);
   }

   // TODO: NOT threadsafe! WriteOutput(0);
}


static int isexcluded(long w)
{
   int i;

   for (i = 0; i < nexclude; ++i) {
      if ((w >= exclude[i][0]) && (w <= exclude[i][1])) {
         return 1;
      }
   }
   return 0;
}


static void parselist(const char *c, long list[][2], int *count)
{
   long a, b;

   while (*c != 0) {
      a = b = 0;
      while (*c >= '0' && *c <= '9') {
         a = a * 10 + (*c++ - '0');
      }
      if (*c == '-') {
         ++c;
         while (*c >= '0' && *c <= '9') {
            b = b * 10 + (*c++ - '0');
         }
      } else {
         b = a;
      }
      if ((a == 0) || (b == 0) || (a > b)) {
         mtxAbort(MTX_HERE,"BAD ARGUMENTS");
      }
      list[*count][0] = a;
      list[*count][1] = b;
      ++*count;
      if (*c == ',') {
         ++c;
      }
   }
}


static void parseCommandLine()
{
   const char *c;

   opt_G = appGetOption(App,"-G --gap");
   opt_n = appGetOption(App,"-n --no-condensation");
   opt_p = appGetOption(App,"-p --use-polynomials");
   opt_t = appGetOption(App,"-t --make-std-basis");
   opt_b = appGetOption(App,"-b --make-ss-basis");
   opt_k = appGetOption(App,"-k --make-pw-kernel");
   while ((c = appGetTextOption(App,"-e --exclude",NULL)) != NULL) {
      parselist(c,exclude,&nexclude);
   }
   while ((c = appGetTextOption(App,"-i --include",NULL)) != NULL) {
      parselist(c,include,&ninclude);
   }
   appGetArguments(App,1,MAX_MODULES);
//   if (opt_G) {
//      MtxMessageLevel = -100;
//   }
}


static void init(int argc, char **argv)
{
   App = appAlloc(&AppInfo,argc,argv);
   parseCommandLine();
   MTX_LOGI("Start pwkond - Peak word condensation");

   #if defined(MTX_DEFAULT_THREADS)
   pthread_mutex_init(&sharedData.mutex, NULL);
   #endif

   loadModules();
   loadConstituents();
   PeakWordsMissing = NumCf;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void addid(Matrix_t *m, FEL f)
{
   int i;
   PTR x;
   if (f == FF_ZERO) {
      return;
   }
   for (i = 0, x = m->data; i < m->nor; ++i, ffStepPtr(&x, m->noc)) {
      ffInsert(x,i,ffAdd(ffExtract(x,i),f));
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////


static void peakWordFound_pex(void* arg)
{
   peakWordFound((struct cf_struct*)arg);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void tryLinear2(long w, FEL f)
{
   int i, ppos = -1;
   long nul;

   for (i = 0; i < NumCf; ++i) { // For each composition factor...
      Matrix_t *word;

      word = wgMakeWord(CfList[i].Wg,w);
      addid(word,f);
      nul = matNullity__(matDup(word));
      if ((nul != 0) && (nul != CfList[i].Info->spl)) {
         matFree(word);
         return;
      }
      if (nul == CfList[i].Info->spl) {
         // possibly a peak word for this constituent
         if (ppos >= 0 || CfList[i].Info->peakWord > 0) {
            matFree(word);
            return;
         }
         nul = matNullity__(matMul(matDup(word),word));
         if (nul != CfList[i].Info->spl) {
            matFree(word);
            return;      // Nullity is not stable
         }
         // This is a peak word candidate for the i-th constituent.
         ppos = i;
      }
      matFree(word);
   }

   if (ppos > -1) { // we have found a new peak word
      struct cf_struct* cf = CfList + ppos;
      cf->Info->peakWord = w;

      // Calculate the nullspace (needed later for standard basei).
      Matrix_t *word = wgMakeWord(cf->Wg, cf->Info->peakWord);
      addid(word, f);
      cf->PWNullSpace = matNullSpace__(word);

      // Compute the peak polynomial (linear case)
      Poly_t *pp = polAlloc(ffOrder,1);
      pp->data[0] = f;
      cf->Info->peakPol = pp;

      --PeakWordsMissing;
      pexExecute(NULL, peakWordFound_pex, cf);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// For a fixed word W, given by its word number, this function finds all peak words of the form
// W+λ1 with λ∈F.

static void tryLinear(long w)
{
   for (uint32_t f = 0; f < ffOrder && PeakWordsMissing > 0; ++f)
   {
      tryLinear2(w,ffFromInt(f));
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int tryp2(long w, int cf, Poly_t *pol)
{
   int i;

   for (i = 0; i < NumCf; ++i) {
      Matrix_t *word, *wordp;
      long nul;

      if (i == cf) {
         continue;
      }
      word = wgMakeWord(CfList[i].Wg,w);
      wordp = matInsert(word,pol);
      matFree(word);
      nul = matNullity__(wordp);
      if (nul != 0) {
         return -1;
      }
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

// For a fixed word W, given by its word number, this function finds all peak words of the form
// p(W) where p∈F[x].

static int tryPoly(long w)
{
   int i;

   for (i = 0; i < NumCf; ++i) {
      Matrix_t *word;
      FPoly_t *mp;

      if (CfList[i].Info->peakWord > 0) {
         continue;                      // We already have a peak word
      }
      word = wgMakeWord(CfList[i].Wg,w);
      mp = minpol(word);
      MTX_LOG2("Constituent %d, minpol = %s", i, fpToEphemeralString(mp));
      int k;
      for (k = 0; k < (int) mp->nFactors; ++k) {
         if (mp->factor[k]->degree * mp->mult[k] == CfList[i].Info->spl) {
            Matrix_t *wp, *wp2;
            long nul;

            MTX_LOG2("%d, factor=%s",i,polToEphemeralString(mp->factor[k]));
            if (tryp2(w,i,mp->factor[k]) == -1) {
               continue;
            }

            // Check if the nullity is stable
            wp = matInsert(word,mp->factor[k]);
            wp2 = matMul(matDup(wp),wp);
            matFree(wp);
            nul = matNullity__(wp2);
            if (nul != CfList[i].Info->spl) {
               continue;
            }
            break;
         }
      }

      if (k < (int) mp->nFactors) {
         CfList[i].Info->peakWord = w;
         CfList[i].Info->peakPol = polDup(mp->factor[k]);
         CfList[i].PWNullSpace = matNullSpace__(matInsert(word,mp->factor[k]));
         --PeakWordsMissing;
         pexExecute(NULL, peakWordFound_pex, CfList + i);
         //peakWordFound(CfList + i);
      } else {
         k = -1;        // Not found
      }
      fpFree(mp);
      matFree(word);
      if (k >= 0) {
         return i;
      }
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void tryWord(long w)
{
   if (isexcluded(w)) {
      return;
   }
   static uint64_t progressTimer = 0;
   if (sysTimeout(&progressTimer, 10))
      MTX_LOGD("Word %ld",w);
   if (opt_p) {
      tryPoly(w);
   } else {
      tryLinear(w);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
   int i;
   long w;

   init(argc, argv);

   if (ninclude > 0) {
      MTX_LOGD("Trying words from inclusion list");
      for (i = 0; PeakWordsMissing > 0 && i < ninclude; ++i) {
         for (w = include[i][0]; w <= include[i][1]; ++w) {
            tryWord(w);
         }
      }
   }
   for (w = 1; PeakWordsMissing > 0; ++w) {
      tryWord(w);
   }

   pexWait();

   WriteOutput(1);
   pexShutdown();
   if (App != NULL) { appFree(App); }
   return 0;
}


/**
@page prog_pwkond pwkond - Peak Word Condensation

@section pwkond_syntax Command Line
<pre>
pwkond @em Options [-Gnptbk] [-i @em List] [-e @em List] @em Name [@em Name ...]
</pre>

@par @em Options
Standard options, see @ref prog_stdopts
@par -G
Produce output in GAP format.
@par -n
Find peak words only, do not condense.
@par -p
Use full polynomials instead of linear factors only in peak word search.
@par -t
Transform generators into standard basis.
@par -b
Calculate a semisimplicity basis.
@par -k
Compute kernel of peak words.
@par -i @em List
Word to try first, for example "-i 100,20-35".
@par -e @em List
Words to be excluded, for example "-e 3,20-99".

@par @em Name
Name of the representation.

@section pwkond_inp Input Files
@par @em Name.cfinfo
Constituent info file.
@par @em Name.1, @em Name.2, ...
Generators.

@par @em NameCF.1, @em NameCF.2, ...
Generators on the irreducible constituents. NameCF is the representation name
followed by the dimension and a letter to distinguish inequivalent constituents
of the same dimension (see @ref prog_chop "chop").

@section pwkond_out Output Files
@par @em Name.cfinfo
Constituent info file.

@par @em NameCF.1k, @em NameCF.2k, ...
Condensed generators on the irreducible constituents.

@par @em NameCF.1.std, @em NameCF.2.std, ...
Condensed generators in standard basis (with -t).

@par @em NameCF.op
Spin-up script for the standard basis (with -t).

@par @em NameCF.np
Condensed peak word.

@par @em NameCF.im
Image used for the condensation.

@par @em NameCF.k
Peak word kernel (with -k or without -n).

@par @em NameCF.ssb
Semisimplicity basis (with -b).


@section pwkond_desc Description


The PWKOND program is part of the @ref sec_progs_lattice "Submodule Lattice Package".
After the irreducible constituents of a module, or a number of modules,
have been found with @ref prog_chop "chop", PWKOND can be used
- to calculate peak words for the constituents,
- to condense the module using the peak words,
- to transform the generators on the constituents to the standard basis as defined
  by the peak word kernel, and
- to calculate a basis reflecting the direct decomposition of the module,
  if the module is semisimple.

By definition, a "peak word" for the i-th constituent is an algebra element W which
fulfills the following conditions:
- W has minimal nullity on the i-th constituent (i.e., its nullity equals the splitting
  field degree for this constituent).
- The nullity is stable, i.e., W and W<sup>2</sup> have the same nullity on the i-th
  constituent.
- W operates regularly (with nullity 0) on all other constituents.

When more than one module is specified on the command line, the peak words found by
@b pwkond are "global", i.e., each peak word selects exactly one of the constituents
of all the modules. Running @b pwkond successively on two modules does not generally
produce global peak words, since a peak word found for module M may have a non-zero
nullity on a different constituent that occurs in another module N but not in M.

The -e option can be used to exclude certain words from the search.
@em List is a list of integers or ranges of integers, for example
"-e 57,82-112,289".
Using "-i" you can specify a list of words which will be tested first.
This can significantly reduce computation time if you already know one
or more peak words for a given module.
The "-n" option disables the condensation phase. If this option is used,
the program stops after the peak words have been found.
If the "-t" option is specified, @b pwkond transforms the generators of all
irreducible constituents to the standard basis defined by the peak word.

For each composition factor there are several output files. If, for
example, one composition factor is X10a, @b pwkond will produce
the following files:
- X10a.std.1 and X10a.std.1 are the operation of the generators on
the constituent with respect to the standard basis defined by the
peak word. These files are created only if the `-t' option is used.
- X10a.op Spin-up script for the standard basis. See ZSB for details.
- X10a.1k and X10a.2k are th action of the generators on the
condensed module.
- X10a.np Condensed peak word. This is a nilpotent matrix.
- X10a.im Image of the peak word.
- X10a.k Kernel of the peak word.
The .cfinfo file is written each time a peak word is found. So, if
the program does not terminate or dies unexpectedly the information about
the peak words found so far are not lost.

If the module is semisimple, @b pwkond can
calculate a basis that respects the decomposition into irreducible
constituents. With respect to this basis, the generators are in block
diagonal form, where the blocks occur in the order determined by @ref prog_chop "chop".
All blocks corresponding to the same constituent are equal, not only
equivalent, and the blocks occur in their "natural" order (as defined by
@ref prog_chop "chop"). This is essential for the tensor condensation procedure
(see @ref prog_precond "precond"). To calculate the semisimplicity basis, use the
"-b" option.
The basis is written to @em Name.ssb. Using "-b" with a module that is not
semisimple produces undefined results. Most probably, @b pwkond will stop
with the error message "row index out of range", or it will write a
singular matrix to @em Name.ssb.

@section pwkond_impl Implementation Details
Internally, a peak word is represented by a pair (n,p) where n is
the canonical number of the word (See @ref prog_zmw "zmw"), and p is a
polynomial. The peak word represented by this pair is p(Wn), Wn
being the n-th word. Without "-p", @b pwkond considers only linear
polynomials. If the "-p" option is used, @b pwkond can find polynomials
of any degree.

Whenever a peak word is found, the generalized condensation
is calculated as follows: The peakword is caculated as a matrix acting on V,
which is then repeatedly raised to higher powers until the nullity stabilizes.
The stable nullity equals the multiplicity k of the constituent times the
degree [E:F] of the splitting field extension.
Having a power w^N of the peakword with stable nullity,
the condensation onto its kernel, i.e., the projection of V onto V/w^N(V),
is determined in the same way as in the @ref prog_zqt "zqt" program.
**/
 
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
