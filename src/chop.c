////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Chop a representation (Calculate composition series)
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <inttypes.h>
#include <string.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constants and macros
////////////////////////////////////////////////////////////////////////////////////////////////////

#define MAX_WORDS 100000    // Maximum number of words to try (for both splitting and idword)
#define MAXFP   6           // Fingerprint size
#define MAXENDO 10          // Max. dimension of endomorphism ring

#define MATFREE(x) { if ((x) != NULL) { matFree(x); (x) = NULL; } }

////////////////////////////////////////////////////////////////////////////////////////////////////
// Types
////////////////////////////////////////////////////////////////////////////////////////////////////

/// A submodule.
typedef struct nodestruct {
   unsigned nodeId;
   char logPrefix[50];          // "[ID:DIM]"
   struct nodestruct *sub, *quot, *parent;
   uint32_t dim;                // Dimension
   int num;                     // (irreducibles only) constituent number
   MatRep_t *Rep;               // Generators
   MatRep_t *TrRep;             // Transposed generators
   long spl;                    // Degree of splitting field
   uint32_t idWord;             // Word used for std basis
   Poly_t *idPol;               // Polynomial "  "    "   
   Poly_t *f1, *f2;             // characteristic polynomial c=f1*f2
   uint32_t fprint[MAXFP];      // Fingerprint
   BitString_t *badWords;
   Matrix_t *nsp;               // Null space 
   uint32_t ggt;                // g.c.d. of nullities
   WgData_t *wg;                // Used by the word generator
   uint32_t wnum;
   Matrix_t *word;
   Charpol_t* cpState;
} node_t;

static void Chop(node_t *n);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////////////////////////////////

static long charpolSeed = 0; 

node_t *root;               // Root node of the constituent tree
long opt_deglimit = -1;     // Max. degree of irred. polynomials
long opt_nullimit;          // Max nullity for exhaustive vector search
node_t *irred[LAT_MAXCF];   // List of irred. constituents
long firstword = 1;
int opt_G = 0;              // GAP output
int opt_i = 0;              // -i: read an existing .cfinfo file
BitString_t *goodWords;     // List of `good' words
static unsigned nodeId = 0; // Counter for node IDs.
LatInfo_t* LI;              // Data for .cfinfo

static long stat_svsplit = 0; // Statistics
static long stat_cpirred = 0;
static long stat_nssplit = 0;
static long stat_dlsplit = 0;
static long stat_irred = 0;
static long stat_exsplit = 0;

static MtxApplicationInfo_t AppInfo = {
   "chop", "Find irreducible constituents",
   "SYNTAX\n"
   "    chop [<Options>] <Name>\n"
   "\n"
   "ARGUMENTS\n"
   "    <Name> .................. Name of the representation\n"
   "\n"
   "OPTIONS\n"
   MTX_COMMON_OPTIONS_DESCRIPTION
   "    -G ...................... GAP output (implies -Q)\n"
   "    -g <NGen> ............... Set number of generators (default is 2)\n"
   "    -s <Word> ............... Start with word number <Word>\n"
   "    -n <MaxNul> ............. Set limit on nullity\n"
   "    -d <MaxDeg> ............. Set limit on degrees of polynomials\n"
   "    -i ...................... Read <Name>.cfinfo, if it exists\n"
   "\n"
   "FILES\n"
   "    <Name>.{1,2,...} ........ I Generators\n"
   "    <Name>.cfinfo ........... O Constituent info file\n"
   "    <Name><Cf>.{1,2...} ..... O Generators on the constituents\n"
};

static MtxApplication_t *App = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////
   
/// Creates a new submodule node.

static node_t *CreateNode(MatRep_t *rep, node_t *parent)
{
   node_t *n = ALLOC(node_t);
   memset(n,0,sizeof(node_t));

   n->nodeId = nodeId++;
   n->parent = parent;
   n->Rep = rep;
   n->dim = rep->Gen[0]->nor;
   snprintf(n->logPrefix, sizeof(n->logPrefix), "[%u:%"PRIu32"]", n->nodeId, n->dim);
   n->num = -1;
   n->spl = -1;
   n->idWord = 0;
   n->badWords = parent ? bsDup(parent->badWords) : bsAllocEmpty();
   n->nsp = NULL;
   n->wg = wgAlloc(n->Rep);
   if (n->wg == NULL) {
      mtxAbort(MTX_HERE,"wgAlloc() failed");
      return NULL;
   }
   n->wnum = -1;
   return n;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Cleans up a submodule node and releases internally allocated memory.
/// 
/// @param n Pointer to the <node_t> structure
/// @param complete Free generators, too. If 0, generators are not freed.

static void cleanUpNode(node_t* n, int complete)
{
   if (n->wg != NULL) { wgFree(n->wg); n->wg = NULL; }
   if (complete) {
      if (n->Rep != NULL) {
         mrFree(n->Rep);
         n->Rep = NULL;
      }
      if (n->idPol != NULL) {
         polFree(n->idPol);
         n->idPol = NULL;
      }
   }
   if (n->badWords != NULL) { bsFree(n->badWords); n->badWords = NULL; }
   if (n->TrRep != NULL) {
      mrFree(n->TrRep);
      n->TrRep = NULL;
   }
   MATFREE(n->nsp);
   MATFREE(n->word);
   if (n->f1 != NULL) { polFree(n->f1); n->f1 = NULL; }
   if (n->f2 != NULL) { polFree(n->f2); n->f2 = NULL; }
   if (n->cpState != NULL) { charpolFree(n->cpState); n->cpState = NULL; }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Makes a word.

static void MakeWord(node_t *n, uint32_t w)
{
   MATFREE(n->word);
   n->word = wgMakeWord(n->wg,w);
   n->wnum = w;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Inserts the current word into a polynomial and calculates the kernel of the resulting matrix.

static void InsertWord(node_t *n, Poly_t *p)
{
   Matrix_t *m = matInsert(n->word,p);
   MATFREE(n->nsp);
   n->nsp = matNullSpace__(m);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Reads the generators and creates the root node of the constituent tree.

static void CreateRoot()
{
   const int scope = mtxBegin(MTX_HERE, "Load module");
   if (opt_i) {
      // Read the number of generators from the cfinfo file.
      LatInfo_t* li = latLoad(LI->BaseName);
      LI->NGen = li->NGen;
      latDestroy(li);
      MTX_LOGD("Set number of generators = %d from existing .cfinfo",LI->NGen);
   }

   MatRep_t *rep = mrLoad(LI->BaseName,LI->NGen);
   LI->field = ffOrder;
   root = CreateNode(rep,NULL);
   mtxEnd(scope);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Returns 1 if the given word is in the "bad words" list for this constituent or any of its
// ancestors, or 0 otherwise.

static int IsBadWord(long w, node_t *n)
{
   if (n == NULL) {
      return 0;
   }
   if (n->badWords != NULL && bsTest(n->badWords,w)) {
      return 1;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Writes the composition series to a string buffer

static void printCompositionSeries(StrBuffer_t* sb, node_t *n, int leading)
{
   MTX_ASSERT(n != NULL);
   if (n->sub == NULL) {        // Irreducible
      int i;
      for (i = 0; LI->Cf[i].dim != n->dim || LI->Cf[i].num != n->num; ++i);
      if (opt_G) {
         printf(leading ? "%d" : ",%d",i + 1);
      } else {
         sbPrintf(sb, " %s", latCfName(LI,i));
      }
   } else {
      printCompositionSeries(sb, n->sub, leading);
      printCompositionSeries(sb, n->quot, 0);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// WriteResult() - Write some information to stdout and create the .cfinfo file.

static void WriteResult(node_t *root)
{
   MTX_LOGI("Chopping completed: %d different composition factors", LI->nCf);
   MTX_LOGI("Writing %s.cfinfo", LI->BaseName);
   latSave(LI);

   // Write composition factors
   if (opt_G) {
      printf("MeatAxe.CompositionFactors := [\n");
      for (int i = 0; i < LI->nCf; ++i) {
         printf("  [ \"%s\", %d, %ld ]",latCfName(LI,i), LI->Cf[i].mult,irred[i]->spl);
         if (i < LI->nCf - 1) {printf(",");}
         printf("\n");
      }
      printf("\n];\n");
   }
   else {
      MTX_LOGI("%s","");
      MTX_LOGI("Name   mult  SF  Fingerprint");
      for (int i = 0; i < LI->nCf; ++i) {
         MTX_XLOGI(sb) {
            sbPrintf(sb, "%-6s %4d  %2ld  ",latCfName(LI,i), LI->Cf[i].mult,irred[i]->spl);
            for (int k = 0; k < MAXFP; ++k) {
               sbPrintf(sb, "%s%lu", k == 0 ? "" : ",", (unsigned long) irred[i]->fprint[k]);
            }
         }
      }
   }

   // Write the composition series
   if (opt_G) {
      StrBuffer_t* sb = sbAlloc(100);
      sbAppend(sb, "MeatAxe.CompositionSeries := [\n");
      printCompositionSeries(sb, root, 1);
      sbAppend(sb, "];\n");
      fputs(sbData(sb),stdout);
      sbFree(sb);

   } else {
      MTX_LOGI("%s","");
      MTX_XLOGI(sb) {
         sbPrintf(sb, "Ascending composition series:");
         printCompositionSeries(sb, root, 1);
      }
   }

   // Write statistics
   MTX_LOGD(" ");
   MTX_LOGD("Statistics:");
   MTX_LOGD("   Saved vectors split: %4ld",stat_svsplit);
   MTX_LOGD("   c(x) irreducible:    %4ld",stat_cpirred);
   MTX_LOGD("   Normal split:        %4ld",stat_nssplit);
   MTX_LOGD("   Dual split:          %4ld",stat_dlsplit);
   MTX_LOGD("   Exceptional split:   %4ld",stat_exsplit);
   MTX_LOGD("   Irreducible:         %4ld",stat_irred);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Splits a constituent
/// 
/// <n>: Pointer to the node to split.
/// <submodule>: Basis for the submodule.
/// <tr>: Indicates that it was a `dual split'
/// 
/// This function is called after a proper submodule has been found. It
/// calculates the action of the generators on both submodule and quotient,
/// and creates two new <note_t> structures for the two parts.
/// The original module, <n>, is cleaned up and <Chop()> is called with
/// both submodule and quotient.

static void splitnode(node_t* n, const Matrix_t* submodule, int tr)
{
   MatRep_t* sub = NULL, * quot = NULL;

   // Split the constituent
   Split(submodule, tr ? n->TrRep : n->Rep, &sub, &quot);

   // If it was a dual split, subspace and quotient have been calculated in the dual module.
   // To get back to the original module, transpose again and exchange sub and quot.
   if (tr) {
      int i;
      Matrix_t* x, * y;
      for (i = 0; i < LI->NGen; ++i) {
         x = matTransposed(sub->Gen[i]);
         matFree(sub->Gen[i]);
         y = matTransposed(quot->Gen[i]);
         matFree(quot->Gen[i]);
         sub->Gen[i] = y;
         quot->Gen[i] = x;
      }
   }

   // Make new nodes for subspace and quotient
   n->sub = CreateNode(sub, n);
   n->quot = CreateNode(quot, n);

   MTX_LOGI("%s Split: Subspace=%u:%"PRIu32", Quotient=%u:%"PRIu32"",
      n->logPrefix, n->sub->nodeId, n->sub->dim, n->quot->nodeId, n->quot->dim);

   // Project saved vectors on the quotient
   if (!tr && n->nsp != NULL) {
      n->quot->nsp = QProjection(submodule, n->nsp);
      matEchelonize(n->quot->nsp); // remove zero vectors
   }

   // Clean up
   cleanUpNode(n, 1);

   // Chop the subspace and quotient
   Chop(n->sub);
   Chop(n->quot);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Finds a vector in 'space' which is not a linear combination of 'basis'.
// 'basis' must be linearly independent.

static Matrix_t* extendbasis(Matrix_t* basis, const Matrix_t* space)
{
   long i, j;
   long dimb = basis->nor;
   long dims = space->nor;
   PTR tmp, x, y;
   FEL f;

   // Concatenate basis and space
   tmp = ffAlloc(dimb + dims, basis->noc);
   memcpy(tmp, basis->data, ffSize(dimb, basis->noc));
   x = (PTR)((char*)tmp + dimb * ffRowSize(basis->noc));
   memcpy(x, space->data, ffSize(dims, basis->noc));

   // Clean with basis
   for (i = 0, x = tmp; i < dimb; ffStepPtr(&x, basis->noc), ++i) {
      const uint32_t piv = ffFindPivot(x, &f, basis->noc);
      if (piv == MTX_NVAL) {
         mtxAbort(MTX_HERE, "extendbasis(): zero vector in basis");
      }
      y = x;
      for (j = i + 1; j < dimb + dims; ++j) {
         ffStepPtr(&y, basis->noc);
         ffAddMulRow(y, x, ffSub(FF_ZERO, ffDiv(ffExtract(y, piv), f)), basis->noc);
      }
   }

   // Find the first non-zero row
   x = ffGetPtr(tmp,dimb,basis->noc);
   for (j = 0; ffFindPivot(x, &f, basis->noc) == MTX_NVAL; ++j, ffStepPtr(&x, basis->noc)) {}
   ffFree(tmp);
   if (j > dims) {
      mtxAbort(MTX_HERE, "extendbasis() failed");
      return NULL;
   }
   return matDupRows(space, j, 1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Checks if a given representation's splitting field has degree [E:F] = dim(V),
/// where V is a given subspace (usually the kernel of an algebra element).
/// Returns 1 if [E:F]=dim(V) or zero otherwise.


static int checkspl(node_t* n, const MatRep_t *rep, Matrix_t *nsp)
{
   int result = 0;

   MTX_LOG2("%s checkspl(): nsp=%lu", n->logPrefix, (unsigned long)nsp->nor);

   // Take the first vector from nsp and change to standard basis.
   Matrix_t* v1 = matDupRows(nsp, 0, 1);
   Matrix_t* const sb1 = spinupStandardBasis(NULL, v1, rep, SF_FIRST);
   MTX_ASSERT(sb1 != NULL && sb1->nor == sb1->noc);
   MatRep_t* const sr1 = mrChangeBasis2(rep, sb1);

   MatRep_t *endo = mrAlloc(0,NULL,0);

   while (1) {

      // Spin up v1 under all endomorphisms found so far. If this yields the whole null-space,
      // we know that the endomorphism ring has at least dimension dim(nsp).
      Matrix_t* subsp = spinup(v1, endo);
      MTX_ASSERT(subsp != NULL);
      if (subsp->nor == nsp->nor) {
         matFree(subsp);
         result = 1;    // Successfull!
         break;
      }

      // Take a vector which is not in «subsp» and make the standard basis again.
      Matrix_t *v2 = extendbasis(subsp,nsp);
      matFree(subsp);
      Matrix_t* sb2 = spinupStandardBasis(NULL, v2,rep, SF_FIRST);
      MTX_ASSERT(sb2 != NULL && sb2->nor == sb2->noc);
      matFree(v2);
      result = mrAreIsomorphic(sr1, rep, sb2);
      if (result != 0)
      {
         // They are identical, i.e., we have found an endomorphism.
         // Put it into the list and try again.
         mrAddGenerator(endo,matMul(matInverse(sb2),sb1),0);
      }
      matFree(sb2);

      if (result == 0) {
         break; // Not successfull
      }
   }

   // Clean up
   matFree(v1);
   matFree(sb1);
   mrFree(sr1);
   mrFree(endo);

   MTX_LOG2("%s checkspl(): result=%d", n->logPrefix, result);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Finds an identifying word ("idWord") for an irreducible constituent.
/// An idWord is an element of the group algebra with minimal nullity, i.e., the nullity equals the
/// splitting field degree [E:F].
/// The word is stored in n->idWord, the polynomial in n->idPol, and its null-space in n->nsp.

static void FindIdWord(node_t* n)
{
   long count = 0;
   const int context = mtxBegin(MTX_HERE, "Searching idword for %s", n->logPrefix);

   // Main loop: Try all words.
   for (uint32_t wordNumber = 1; n->idWord == 0; ++wordNumber) {
      static uint64_t progressTimer = 0;
      if (sysTimeout(&progressTimer, 5)) {
         MTX_LOGD("%s Searching idWord (%"PRIu32")...", n->logPrefix, wordNumber);
      }
      if (IsBadWord(wordNumber, n)) {
         continue;
      }

      // Make the word and its characteristic polynomial
      MakeWord(n, wordNumber);
      MTX_LOG2("%s Next word: %"PRIu32", gcd=%"PRIu32, n->logPrefix, wordNumber, n->ggt);
      FPoly_t* cpol = charpol(n->word);
      if (charpolSeed >= n->word->nor) { charpolSeed = 0; } // TODO: remove (changes tests)
      MTX_LOG2("%s c(x)=%s", n->logPrefix, fpToEphemeralString(cpol));
      for (uint32_t k = 0; k < cpol->nFactors; ++k) {
         n->ggt = gcd(n->ggt, cpol->mult[k] * cpol->factor[k]->degree);
      }

      // Try all factors of c(x) with degree <= g.c.d. of all degrees
      for (uint32_t k = 0; k < cpol->nFactors; ++k) {
         if (cpol->factor[k]->degree > n->ggt)
            continue;
         if (++count > MAX_WORDS)
            mtxAbort(MTX_HERE, "FindIdWord() failed");
         InsertWord(n, cpol->factor[k]);
         n->ggt = gcd(n->ggt, n->nsp->nor);
         if (n->nsp->nor > n->ggt)
            continue;
         MTX_XLOG2(msg) {
            sbPrintf(msg, "%s factor=", n->logPrefix);
            polFormat(msg, cpol->factor[k]);
            sbPrintf(msg, ", nsp=%"PRIu32", gcd=%"PRIu32, n->nsp->nor, n->ggt);
         }

         if (checkspl(n, n->Rep, n->nsp)) {
            bsSet(goodWords, wordNumber);
            n->idWord = wordNumber;
            n->idPol = polDup(cpol->factor[k]);
            break;
         }
      }
      fpFree(cpol);
   }

   MTX_XLOGD(msg) {
      sbPrintf(msg, "%s idWord=%"PRIu32", idPol=", n->logPrefix, n->idWord);
      polFormat(msg, n->idPol);
   }
   mtxEnd(context);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Checks if a given irreducible module is already contained in the list of composition factors.
// If yes, returns the index. If not, inserts the new irreducible and return its index.

static void newirred(node_t* n)
{
   int i;
   Matrix_t* b;

   // Check if the module is already in the list
   wgMakeFingerPrint(n->wg, n->fprint);
   for (i = 0; i < LI->nCf && n->dim >= irred[i]->dim; ++i) {
      // Compare dimensions and fingerprints
      if (n->dim != irred[i]->dim ||
          memcmp(n->fprint, irred[i]->fprint, sizeof(n->fprint))) {
         continue;
      }

      if (IsIsomorphic(irred[i]->Rep, LI->Cf + i, n->Rep, NULL, 0)) {
         ///++irred[i]->mult;
         ++LI->Cf[i].mult;
         n->num = irred[i]->num;
         MTX_LOGI("%s Irreducible (%s)", n->logPrefix, latCfName(LI, i));
         cleanUpNode(n, 1);
         return;
      }
   }

   // It's a new irreducible!
   if (LI->nCf >= LAT_MAXCF) {
      mtxAbort(MTX_HERE, "TOO MANY CONSTITUENTS");
   }
   for (int k = LI->nCf - 1; k >= i; --k) {
      irred[k + 1] = irred[k];
      LI->Cf[k + 1] = LI->Cf[k];
   }
   irred[i] = n;
   ++LI->nCf;

   // Set the fields
   if (i == 0 || irred[i]->dim != irred[i - 1]->dim) {
      irred[i]->num = 0;
   }
   else {
      irred[i]->num = irred[i - 1]->num + 1;
   }
   LI->Cf[i].mult = 1;

   // Make idWord and change to std basis
   MATFREE(n->nsp);
   MTX_ASSERT(n->idWord == 0);
   FindIdWord(n);
   LI->Cf[i].dim = irred[i]->dim;   // required for cfname()
   LI->Cf[i].num = irred[i]->num;
   LI->Cf[i].mult = 1;
   LI->Cf[i].idWord = n->idWord;
   LI->Cf[i].idPol = polDup(n->idPol);
   LI->Cf[i].spl = n->spl = n->nsp->nor;
   b = spinupStandardBasis(NULL, n->nsp, n->Rep, SF_FIRST);
   MTX_ASSERT(b != NULL && b->nor == b->noc);
   mrChangeBasis(n->Rep, b);
   matFree(b);
   MTX_LOGI("%s Irreducible (%s)", n->logPrefix, latCfName(LI, i));

   // Write out the generators
   if (irred[i]->spl > 1) {
      MTX_LOGD("%s Splitting field has degree %ld", n->logPrefix, irred[i]->spl);
   }
   for (int k = 0; k < LI->NGen; ++k) {
      char fn[200];
      sprintf(fn, "%s%s.%d", LI->BaseName, latCfName(LI, i), k + 1);
      matSave(irred[i]->Rep->Gen[k], fn);
   }
   cleanUpNode(n, 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
   
// Tries to split a module using saved vectors.
//   Returns 1 on success, 0 otherwise.

static int SplitWithSavedVectors(node_t* n)
{
   if (n->nsp == NULL || n->nsp->nor == 0) {
      return 0;
   }
   MTX_LOGD("%s Trying %lu saved vectors", n->logPrefix, (unsigned long) n->nsp->nor);

   Matrix_t* submodule = NULL;
   if (spinupFindSubmodule(&submodule, n->nsp, n->Rep, SF_EACH, 0) > 0) {
      MTX_LOGD("%s Splitting with saved vectors succeeded", n->logPrefix);
      ++stat_svsplit;
      splitnode(n, submodule, 0);
      matFree(submodule);
      return 1;
   }
   
    MTX_LOGD("%s Splitting with saved vectors failed", n->logPrefix);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Calculates the image of a single vector, or a set of vectors, under p(A).

Matrix_t *polymap(Matrix_t *v, Matrix_t *m, Poly_t *p)
{
   Matrix_t *result, *tmp;
   int i;

   result = matAlloc(v->field,v->nor,v->noc);
   tmp = matDup(v);
   for (i = 0; i <= p->degree; ++i) {
      FEL f = p->data[i];
      PTR x = result->data, y = tmp->data;
      long i;
      for (i = v->nor; i > 0; --i) {
         ffAddMulRow(x,y,f,v->noc);
         ffStepPtr(&x,v->noc);
         ffStepPtr(&y,v->noc);
      }
      matMul(tmp,m);
   }
   matFree(tmp);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
   
// Calculates a vector in the null-space of p(A). p(x) is assumed to be a factor of f1(x), i.e.,
// p(x) must occur in the first cyclic subspace. 

static void make_kern(node_t *n, Poly_t *p)
{
   Poly_t *cof, *f;
   Matrix_t *seed;

   f = polDup(n->f1);
   cof = polDivMod(f,p);
   MTX_ASSERT(f->degree == -1);
   polFree(f);

   if (n->nsp != NULL) {
      matFree(n->nsp);
   }
   seed = matAlloc(ffOrder,1,n->dim);
   ffInsert(seed->data,charpolSeed,FF_ONE);
   n->nsp = polymap(seed,n->word,cof);
   matFree(seed);
   polFree(cof);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Tries to find a vector in the null-space of p(A^T).
/// Returns the vector or NULL if not found.
/// We consider only the first cyclic subspace, so make_trkern() may not find a vector.

static Matrix_t* make_trkern(node_t *n, Poly_t *p)
{
   Matrix_t* result = NULL;

   Matrix_t* mt = matTransposed(n->word);
   Charpol_t* state = charpolStart(mt, PM_CHARPOL, charpolSeed);
   Poly_t* pt = charpolFactor(state); // factor of c(x)
   Poly_t* cofactor = polDivMod(pt,p);
   if (pt->degree == -1) {
      // p divides pt
      Matrix_t *seed = matAlloc(ffOrder,1,n->dim);
      ffInsert(seed->data,charpolSeed,FF_ONE);
      result = polymap(seed,mt,cofactor);
      matFree(seed);
   }
   polFree(cofactor);
   polFree(pt);
   charpolFree(state);
   matFree(mt);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Calculates the characteristic polynomial on the first cyclic subspace.

static FPoly_t *make_f1(node_t *n)
{
   FPoly_t *cpol;

   if (n->f1 != NULL) { polFree(n->f1); n->f1 = NULL;}
   if (n->f2 != NULL) { polFree(n->f2); n->f2 = NULL;}

   if (n->cpState) charpolFree(n->cpState);
   n->cpState = charpolStart(n->word, PM_CHARPOL, charpolSeed);

   n->f1 = charpolFactor(n->cpState);
   if (charpolSeed >= n->word->nor) charpolSeed = 0;    // TODO: remove (changes tests?)

   cpol = Factorization(n->f1);

   // Sort the factors by ascending multiplicity.
   for (uint32_t i = 0; i < cpol->nFactors; ++i) {
      for (uint32_t k = i + 1; k < cpol->nFactors; ++k) {
         long val1 = cpol->mult[i];
         long val2 = cpol->mult[k];
         if (val1 > val2) {
            Poly_t *p = cpol->factor[i];
            long e = cpol->mult[i];
            cpol->factor[i] = cpol->factor[k];
            cpol->mult[i] = cpol->mult[k];
            cpol->factor[k] = p;
            cpol->mult[k] = e;
         }
      }
   }
   MTX_XLOG2(msg) {
      sbPrintf(msg, "[%u:%"PRIu32"] f1(x) = ", n->nodeId, n->dim);
      fpFormat(msg, cpol);
   }
   return cpol;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Complete the characteristic polynomial.
/// This function completes the calculation of the characteristic polynomial. We assume that the
/// first factor (i.e. the characteristic polynomial on the first cyclic subspace) has already been
/// calculated. The remaining factors are stored in «n->f2».

static void make_f2(node_t* n)
{
   if (n->f2 != NULL) {
      return;
   }
   n->f2 = polAlloc(n->f1->field, 0);
   Poly_t* f;
   while ((f = charpolFactor(n->cpState)) != NULL) {
      polMul(n->f2, f);
      polFree(f);
   }
   MTX_XLOG2(msg) {
      sbPrintf(msg, "%s f2(x) = ", n->logPrefix);
      FPoly_t* x = Factorization(n->f2);
      fpFormat(msg, x);
      fpFree(x);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Tries to split the module with the first vector in «n->nsp».
/// Returns 1 on success, 0 on failure

static int SplitWithNsp(node_t* n)
{
   MTX_LOG2("%s Trying to split with null-space", n->logPrefix);

   Matrix_t* v1 = matDupRows(n->nsp, 0, 1);
   Matrix_t* sub = spinup(v1, n->Rep);
   matFree(v1);

   const int haveSubmodule = sub->nor > 0 && sub->nor < sub->noc;
   if (haveSubmodule) {
      ++stat_nssplit;
      bsSet(goodWords, n->wnum);
      splitnode(n, sub, 0);
   } else {
      MTX_LOG2("%s Failed", n->logPrefix);
   }
   matFree(sub);
   return haveSubmodule;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int polMultiplicity(const Poly_t *factor, const Poly_t *pol)
{
   int mult = 0;
   Poly_t *f2 = polDup(pol);
   int done = 0;
   while (!done) {
      Poly_t *q = polDivMod(f2,factor);
      if (f2->degree == -1) {
         ++mult;
         polFree(f2);
         f2 = q;
      } else {
         polFree(q);
         done = 1;
      }
   }
   polFree(f2);
   return mult;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Tries to split a submodule with one irreducible factor of f1(x), where c(x)=f1(x)*f2(x) is the
/// characteristic polynomial of the current word.
///
/// @param pol The irreducible factor of c(x)
/// @param multiplicityInF1 Multiplicity of p in f1(x)
///
/// @return: 1 = success, 0 = failed

static int try_poly(node_t* n, Poly_t* pol, long multiplicityInF1)
{
   Matrix_t* sub;
   int mult;

   // Try to split with a single vector in the null space of p(A)
   make_kern(n, pol);
   if (SplitWithNsp(n)) {
      return 1;
   }

   // Now find out if we can prove irreducibility with Norton's criterion. We need that the
   // factor we just tried occurs with multiplicity 1 in the characteristic polynomial, c(x).
   int canProveIrreducibility = 1;
   if (n->f2 == NULL) {
      make_f2(n);
   }
   mult = polMultiplicity(pol, n->f2) + multiplicityInF1;     // Multiplicity in c(x)

   // If the multiplicity is not one, we can still prove irreducibility, but we must spin up every
   // vector in the null-space of p(A). We do this only if the null-space is small (option -n).
   if (mult > 1) {
      if (mult * pol->degree > opt_nullimit) {
         canProveIrreducibility = 0;
         MTX_LOG2("%s Cannot check for irreducibility, null-space=%d",
            n->logPrefix,
            mult * pol->degree);
      }
      else {
         matFree(n->nsp);
         n->nsp = matNullSpace__(matInsert(n->word, pol));
         MTX_LOG2("%s 2nd spin-up, null-space = %lu", n->logPrefix, (unsigned long)n->nsp->nor);

         Matrix_t* sub2 = NULL;
         if (spinupFindSubmodule(&sub2, n->nsp, n->Rep, SF_MAKE, 0) > 0) {
            ++stat_nssplit;
            bsSet(goodWords, n->wnum);
            splitnode(n, sub2, 0);
            matFree(sub2);
            return 1;
         }
      }
   }

   // Not split, try dual split
   Matrix_t* trNullVector = make_trkern(n, pol);
   if (trNullVector == NULL) {
      MTX_LOG2("%s No seed vector found, dual split skipped", n->logPrefix);
      return 0;
   }
   MTX_LOG2("%s Try dual split...", n->logPrefix);
   if (n->TrRep == NULL) {
      n->TrRep = mrTransposed(n->Rep);
   }
   sub = spinup(trNullVector, n->TrRep);
   matFree(trNullVector);
   const int haveSubmodule = sub->nor > 0 && sub->nor < sub->noc;
   MTX_LOG2("%s Dual split %s", n->logPrefix, haveSubmodule ? "successful" : "failed");
   if (haveSubmodule) {
      ++stat_dlsplit;
      splitnode(n, sub, 1);
   }
   matFree(sub);
   if (haveSubmodule) {
      return 1;
   }

   if (canProveIrreducibility) {
      // The module is irreducible
      newirred(n);
      bsSet(goodWords, n->wnum);
      ++stat_irred;
      return 1;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Try to split in exceptional cases. Returns 1 on success (split), 0 otherwise.
// @p n is a pointer to the module,
// @p cp is the characteristic polynomial of the current word,
// @p cpf is the factorized form of the characteristic polynomial, and
// @p i is the index of the factor to try
//
// This function is called from try_exceptional() for each irreducible
// factor of the characteristic polynomial of degree >= 2.

static int try_ex_factor(node_t* n, Poly_t* cp, FPoly_t* cpf, int factor)
{
   int k;
   long rndword;
   Poly_t* p, * q, * tmp, * gcd[3];
   Matrix_t* B, * iA, * v, * v1;

   MTX_XLOG2(msg) {
      sbAppend(msg, "Trying factor (");
      polFormat(msg, cpf->factor[factor]);
      sbPrintf(msg, ")^%lu", (unsigned long) cpf->mult[factor]);
   }

   // Calculate p(x) = maximal power of the irred factor in c(x)
   p = polDup(cpf->factor[factor]);
   for (k = 1; k < cpf->mult[factor]; ++k) {
      polMul(p, cpf->factor[factor]);
   }

   // Calculate the complement q(x) with q(x)p(x)=c(x)
   tmp = polDup(cp);
   q = polDivMod(tmp, p);
   MTX_ASSERT(tmp->degree == -1);
   polFree(tmp);

   // Calculate i(x):=b(x)q(x) with a(x)p(x)+b(x)q(x) = 1.
   polGcdEx(p, q, gcd);
   MTX_ASSERT(gcd[0]->degree == 0);
   polMul(q, gcd[2]);

   // Insert the word into i(x) and clean up polynomials.
   iA = matInsert(n->word, q);
   polFree(gcd[0]);
   polFree(gcd[1]);
   polFree(gcd[2]);
   polFree(p);
   polFree(q);

   // Choose a second random word, B, and calculate [A,i(A)Bi(A)]
   rndword = n->wnum + mtxRandomInt(42);
   MTX_LOG2("%s Choosing random word %ld", n->logPrefix, rndword);
   B = wgMakeWord(n->wg, rndword);

   // Select a random vector in the image of the commutator
   v = matAlloc(B->field, 1, B->noc);
   for (k = 0; k < v->noc; ++k) {
      ffInsert(v->data, k, ffFromInt(mtxRandomInt(ffOrder)));
   }

   v1 = matDup(v);
   matMul(matMul(matMul(matMul(v, n->word), iA), B), iA);
   matMul(matMul(matMul(matMul(v1, iA), B), iA), B);
   matAddMul(v, v1, ffNeg(FF_ONE));
   matFree(v1);
   matFree(iA);
   matFree(B);

   // Try to split with this vector
   // OLD: Matrix_t* sub = SpinUp(v, n->Rep, SF_FIRST | SF_SUB, NULL, NULL);
   Matrix_t* sub = spinup(v, n->Rep);
   matFree(v);

   const int haveSubmodule = sub->nor > 0 && sub->nor < sub->noc;
   MTX_LOG2("%s Split (exceptional): %s", n->logPrefix, haveSubmodule ? "successful" : "failed");
   if (haveSubmodule) {
      ++stat_exsplit;
      splitnode(n, sub, 0);
   }
   matFree(sub);
   return haveSubmodule;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Try to split in exceptional cases
// Returns 1 on success (module was split), 0 otherwise.
//
// This function tries to split a module. It uses an algorithm developed by G. Ivanyos and K. Lux
// which is specifically designed for the `exceptional' cases where the standard methods fail.

static int try_exceptional(node_t *n)
{
   Poly_t *cp;
   FPoly_t *cpf;
   int success = 0;
   MTX_ASSERT(n->f1 != NULL && n->f2 != NULL);

   MTX_LOG2("%s Trying exceptional cases", n->logPrefix);

   // Calculate the complete characteristic polynomial c(x) and its irreducible factors.
   cp = polDup(n->f1);
   polMul(cp,n->f2);
   cpf = Factorization(cp);

   // Try all factors of degree >= 2
   for (uint32_t factor = 0; factor < cpf->nFactors && !success; ++factor) {
      if (cpf->factor[factor]->degree < 2) {
         continue;
      }
      success = try_ex_factor(n,cp,cpf,factor);
   }

   fpFree(cpf);
   polFree(cp);
   return success;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Tries to chop a module using a given word. Returns 1 on success, 0 otherwise.

static int ChopWithWord(node_t *n, long wn, int try_ex)
{
   long dlimit = opt_deglimit;          // Limit on degree

   MTX_LOG2("%s Next word is %ld (=%s)",n->logPrefix, wn,wgSymbolicName(n->wg,wn));
   charpolSeed = (charpolSeed + 2) % n->dim;    // BUG: '+2' for 2.3 compatibility
   MakeWord(n,wn);
   FPoly_t* f1 = make_f1(n);           // Make first part of c(x)

   // If c(x) is irreducible, then the module is irreducible
   if (f1->factor[0]->degree == n->dim) {
      MTX_LOG2("%s c(x) is irreducible", n->logPrefix);
      newirred(n);
      ++stat_cpirred;
      bsSet(goodWords,wn);
      fpFree(f1);
      return 1;
   }

   // Try all factors of f1(x)
   uint32_t pi;
   int done = 0;
   for (pi = 0; !done && pi < f1->nFactors; ++pi) {
      MTX_XLOG2(msg) {
         sbPrintf(msg, "%s Next factor: (", n->logPrefix);
         polFormat(msg, f1->factor[pi]);
         sbPrintf(msg, ")^%lu", (unsigned long) f1->mult[pi]);
      }
      if (dlimit > 0 && f1->factor[pi]->degree > dlimit) {
         MTX_LOG2("%s deg > %ld -- discarded", n->logPrefix, dlimit);
         continue;
      }
      done = try_poly(n,f1->factor[pi],f1->mult[pi]);
      MTX_LOG2("%s try_poly()=%d", n->logPrefix, done);
   }
   fpFree(f1);
   if (!done && try_ex) {
      done = try_exceptional(n);
   }
   return done;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Handle 1-dimensional modules.
/// If the dimension is 1, marks the constituent as irrreducible and returns 1. 
/// Returns 0 otherwise.

static int IsOneDimensional(node_t *n)
{
   if (n->dim > 1) {
      return 0;
   }
   MTX_LOGD("%s Dimension is one -- irreducible", n->logPrefix);
   newirred(n);
   ++stat_irred;
   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int tryWord(node_t* n, size_t wordNo, int try_ext)
{
   if (IsBadWord(wordNo, n)) {
      MTX_LOG2("%s Skip bad word %ld", n->logPrefix, wordNo);
      return 0;
   }
   MTX_LOG2("%s Trying word %ld", n->logPrefix, wordNo);
   if (ChopWithWord(n, wordNo, try_ext)) {
      return 1;
   }
   MTX_LOG2("%s Add bad word %ld", n->logPrefix, wordNo);
   bsSet(n->badWords, wordNo);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Chops a constituent or proves that it is irreducible.

static void Chop(node_t *n)
{
   int count = 0;

   if (n == NULL) {
      mtxAbort(MTX_HERE,"node=NULL: %s",MTX_ERR_BADARG);
   }
   MTX_LOGI("%s Chop: Dim=%"PRIu32"", n->logPrefix, n->dim);

   // Handle dimension 1.
   if (IsOneDimensional(n)) {
      return;
   }
   // Try splitting with saved vectors.
   if (SplitWithSavedVectors(n)) {
      return;
   }

   MTX_LOG2("%s Trying known good words", n->logPrefix);
   {
      size_t wordNo = 0;
      if (bsFirst(goodWords, &wordNo)) {
         do {
            if (tryWord(n, wordNo, count > 10))
               return;
            ++count;
         } while (bsNext(goodWords, &wordNo));
      }
   }

   MTX_LOG2("%s Trying other words", n->logPrefix);
   for (long wordNo = firstword; count < MAX_WORDS; ++count, ++wordNo) {
      
      if (bsTest(goodWords, wordNo)) {
         continue;
      }
      if (tryWord(n, wordNo, count > 10))
         return;
      bsSet(n->badWords, wordNo);
   }
   
   mtxAbort(MTX_HERE, "GAME OVER");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
   App = appAlloc(&AppInfo,argc,argv);
   const int scope = mtxBegin(MTX_HERE, "Initialize program");
   goodWords = bsAllocEmpty();
   opt_G = appGetOption(App,"-G --gap");
   opt_i = appGetOption(App,"-i --read-cfinfo");
   firstword = appGetIntOption(App,"-s", 1, 1, 100000);
   int ngen = appGetIntOption(App,"-g --generators",2,1,MAXGEN);
   opt_deglimit = appGetIntOption(App,"-d --max-polynomial-degree", -1, -1, 100);
   opt_nullimit = appGetIntOption(App,"-n", 3, 1, 20);
   appGetArguments(App,1,1);
   LI = latCreate(App->argV[0]);
   LI->NGen = ngen;
   mtxEnd(scope);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void cleanup()
{
   bsFree(goodWords);
   for (size_t i = 0; i < LI->nCf; ++i) {
      cleanUpNode(irred[i], 1);
      irred[i] = NULL;
   }
   latDestroy(LI);

   if (App != NULL) {
      appFree(App);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   init(argc,argv);
   MTX_LOGI("Start chop - Find irredicible constituents");
   CreateRoot();
   Chop(root);
   WriteResult(root);
   cleanup();
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT_OFF*

/**
@page prog_chop chop - Find Irreducible Constituents

@section chop_syntax Command Line
<pre>
chop [@em Options] [-Gi] [-g @em NGen] [-s @em Word] [-n @em MaxNul] [-d @em MaxDeg] @em Name
</pre>

@par @em Options
Standard options, see @ref prog_stdopts

@par @em -G --gap
Produce GAP output. Implies -Q.

@par @em -i --read-cfinfo
Read @em Name.cfinfo, if it exists, to determine the number of generators.

@par -g @em NGen
Set the number of generators. Default is two generators, but see -i.

@par -s @em Word
Start with the given word number instead of 1.

@par -n @em MaxNul
Set limit on nullity. Only null-spaces with a dimension less than or equal to
@em MaxNul are searched completely.

@par -d @em MaxDeg
Set limit on degrees of polynomials.

@par @em Name
Name of the module to chop.

@section chop_inp Input Files
@par @em Name.1, @em Name.2, ...
Action of the generators on the module.

@section chop_out Output Files

@par @em Name.cfinfo
Constituent information file

@par @em Name.X.1, @em Name.X.2, ...
Action of the generators on the constituent X.

@section chop_desc Description
The CHOP program is part of the @ref sec_progs_lattice "Submodule Lattice Package".
CHOP calculates the irreducible constituents of a given matrix representation.
The representing matrices of the generators are read from input files, see "input Files"
above. Unless a different number of generators has been specified with -g, two
generators are expected. However, if the "-i" option is used, and the file @em Name.cfinfo
exists, @b chop takes the number of generators from this file and ignores the "-g" option.

For each composition factor @b chop writes the action of the generators to @em CFName.1,
@em CFName.2, ... @em CFName is the name of the composition factor, which is constructed
by appending the dimension and a letter to the module name. For example, "X10a.1"
is the action of the first generator on the first composition factor of dimension
10 of the module X. If a second, inequivalent composition factor of dimension 10
was found, it would be named `X10b' and so on.
@b chop also creates the file @em Name.cfinfo' containing a list of all composition factors.
This file is used by subsequent programs such as @ref prog_pwkond "pwkond".

@section chop_impl Implementation Details
@b chop repeatedly splits a module into submodule and quotient until it arrives at the
irreducible constituents. Thus, it finds a composition series. The program assumes that the
algebra generated by the input matrices contains the unit matrix.

In order to split a given module or to prove its irreducibility the algorithm needs an element
of the algebra with a non-trivial but low-dimensional kernel. Such elements are searched by
taking linear combinations of certain products of the generators ("words").
See the description of the @ref prog_zmw "zmw" program for more details on the word generator.
By default, @b chop tries all words in the order defined by the word generator. The "-s" option
may be used to make @b chop start with a word different from 1.

For each word A generated in this way, the program calculates its characteristic polynomial and
examines the irreducible factors.  If p(x) is an irreducible factor, p(A) has a non-trivial
kernel.  Then, one vector of the kernel is chosen and the smallest submodule containing this
vector is calculated.
If the vector spans a proper submodule, the action of the generators on this submodule as well
as on the quotient are calculated and the same procedure is applied recursively to both
submodule and quotient.

To avoid expensive matrix multiplications in the calculation of p(A),
there is a limit on the degree of p(x). This limit can be set with
the "-d" option and defaults to 5.

If a module cannot be split by the program, it may be irreducible.
In order to prove this, @b chop uses Norton's criterion. This requires,
however, to find an algebra element with a small kernel, because up
to scalar multiples each vector in the kernel must be examined to see
whether it spins up to the whole module. For this reason a "nullity
threshold" m is maintained by the program. Initially, m is set to 3 or
to the value given in the "-n" option.
Each algebra element that has a nullity less then or
equal to m is used for the Norton test.

In some cases the algorithm described is not able to split the module
although it is reducible. These exceptional cases are treated with
an alternative strategy described in @ref LI98 "[LI98]".

Algebra elements with trivial kernel are useless for the algorithm, so
an attempt is made to avoid unnecessary computation of such elements.
Once an element is known to have a trivial kernel on a given module
M, the program will mark it as invertible and ignore it for all
constituents of M.

If a constituent is irreducible but not absolutely irreducible, the
nullity of any element in the algebra will be a multiple of [E:F],
where F is the ground field and E the splitting field.
This situation is recognized by calculating the greatest common divisor d
of all nullities which occur during the search.
In order to prove that the splitting field degree is equal to
@em d, the following method is used: Take a word with nullity @em d
and two vectors v1, v2 in its null-space. Use these vectors as
seeds for a standard basis algorithm. If the resulting representations
are different, [E:F] is less than @em d, and the word is discarded.
Otherwise, the linear map
which transforms one standard basis into the other is an endomorphism
@em e of the module. If v1, under the action of e, spins up to
the whole null space, then [E:F]=@em d. Otherwise, take a third vector
not in the span and repeat the procedure above. Again, this yields an
endomorphism, or it turns out that [E:F]<@em d.
These steps are repeated until a word with nullity [E:F] is found.
 */

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
