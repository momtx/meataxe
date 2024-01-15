////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Spin-up
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <inttypes.h>
#include <string.h>

/// @defgroup g_spinup Spin-Up and Split
/// @{
/// @details
/// Given a matrix representation and a seed vector v, the spin-up algorithm
/// calculates the submodule generated by the seed vector, i.e., the smallest 
/// subspace containing v which is invariant under the generators.
/// SpinUp() can handle multiple seed vectors, search for cyclic vectors
/// generating the whole space, and generate seed vectors as linear combinations
/// of a given basis.
///
/// @anchor spinupscripts <b>Spin-up Scripts</b><br>
/// When spinning up a seed vector, you can record the operations performed
/// by the algorithm in a spin-up script. This script can then be fed into
/// spinupWithScript() to repeat the procedure with a different seed vector 
/// and different generators.
/// A spin-up consists of a pair (v,g)∈ℤ<sup>2</sup> for each subspace basis vector that was
/// produced by the spin-up process. 
/// The pair is either of the form (v,-1) with v>0, meaning that the basis vector is the
/// v-th seed vector. Otherwise the pair is of the form (v,g) with g≥0, meaning that the
/// vector is the image of the v-th basis vector under the g-th generator. Note that
/// vector numbers are 1-based, but generator indices start with 0!
///
/// @anchor stdbasis <b>Standard Basis</b><br>
/// Normally, the basis vectors computed during the spin-up process are chosen randomly.
/// However, the spin-up algorithm can be used in "standard basis" mode.
/// In this mode, the result is invariant under a change of basis.
/// More precisely, if a given seed vector v and generators g<sub>1</sub>,...g<sub>n</sub>
/// produce the basis (b<sub>1</sub>,...b_<sub>m</sub>), and A is a nonsingular matrix,
/// then vA and A<sup>-1</sup>g<sub>1</sub>A,...A<sup>-1</sup>g<sub>n</sub>A produce the basis
/// (b<sub>1</sub>A,...b<sub>m</sub>A).
/// The standard basis algorithm performs a normal spinup but maintains a second
/// basis by mirroring all operations on the original vectors except the cleaning passes.
/// It consumes more memory and computing time, and the resulting basis is not in
/// echelon form.
///
/// <b>Splitting a representation</b><br>
/// If a proper invariant subspace U<V has been found for a matrix representation M,
/// the restriction of M to U as well as the representation on V/U can be calculated.
/// This is called splitting the representation. The basis for V/U is chosen randomly.

#if defined(MTX_DEFAULT_THREADS)
   #define MUTEX_INIT(mutex) pthread_mutex_init(&mutex, NULL)
   #define MUTEX_DESTROY(mutex) pthread_mutex_destroy(&mutex)
   #define MUTEX_LOCK(mutex) pthread_mutex_lock(&mutex)
   #define MUTEX_UNLOCK(mutex) pthread_mutex_unlock(&mutex)
#else
   #define MUTEX_INIT(mutx)
   #define MUTEX_DESTROY(mutx)
   #define MUTEX_LOCK(mutex)
   #define MUTEX_UNLOCK(ctx)
#endif

/// @private
struct Workspace {
    struct Workspace* next;
    int isInPool;
    int isValid;
    PTR rows;
    PTR stdRows;        // Only used with SF_STD, NULL otherwise
    uint32_t* piv;
    int32_t *ops;       // Only used if reccording a spin-up script, NULL otherwise
    uint32_t field;
    uint32_t dim;       // Subspace dimension (= number of valid entries in rows and piv)
    uint32_t noc;       // Row size (number of columns, maximum subspace dimension)
};

/// @private
typedef struct SpinupContext {

   // Constant during context lifetime. Pointers are not owning.
   unsigned flags;
   uint32_t field;		// Field order
   uint32_t noc;		// Dimension of the whole space
   const Matrix_t* seed;
   const MatRep_t* rep;
   uint32_t maxSubspaceDim;

   #if defined(MTX_DEFAULT_THREADS)
   pthread_mutex_t mutex;
   #endif

   Matrix_t* submodule;         // Submodule found by spinup
   IntMatrix_t* script;         // Spin-up script (if requested)
   uint32_t resultSeed;         // Seed vector which produced the result

   struct Workspace* workspaces;// Workspace pool
   PexGroup_t* pexGroup;
} SpinupContext_t;

////////////////////////////////////////////////////////////////////////////////////////////////////

void imatMove(IntMatrix_t** dst, IntMatrix_t** src)
{
   MTX_ASSERT(src != NULL);
   if (dst != NULL) {
      if (*dst != NULL)
         imatFree(*dst);
      *dst = *src;
      *src = NULL;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void matMove(Matrix_t** dst, Matrix_t** src)
{
   MTX_ASSERT(src != NULL);
   if (dst != NULL) {
      if (*dst != NULL)
         matFree(*dst);
      *dst = *src;
      *src = NULL;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static struct Workspace* provideWorkspace(struct SpinupContext* ctx)
{
   struct Workspace* ws = NULL;
   MUTEX_LOCK(ctx->mutex);
   if (ctx->workspaces != NULL) {
      ws = ctx->workspaces;
      ctx->workspaces = ws->next;
   }
   MUTEX_UNLOCK(ctx->mutex);
   if (ws == NULL) {
      ws = ALLOC(struct Workspace);
      ws->rows = ffAlloc(ctx->noc + 1, ctx->noc);
   }
   if (ctx->flags & SF_STD) {
      ws->stdRows = ffAlloc(ctx->noc + 1, ctx->noc);
      ws->ops = NALLOC(int32_t, 2 * (ctx->noc + 1));
   }
   ws->piv = NALLOC(uint32_t, ctx->noc + 1);

   ws->isInPool = 0;
   ws->next = NULL;
   ws->dim = 0;
   ws->noc = ctx->noc;
   ws->field = ctx->field;
   return ws;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void destroyWorkspace(struct Workspace* ws)
{
   MTX_ASSERT(!ws->isInPool);
   sysFree(ws->ops);
   ws->ops = NULL;
   sysFree(ws->piv);
   ws->piv = NULL;
   if (ws->rows != NULL) {
      ffFree(ws->rows);
      ws->rows = NULL;
   }
   if (ws->stdRows != NULL) {
      ffFree(ws->stdRows);
      ws->stdRows = NULL;
   }
   sysFree(ws);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void returnWorkspaceToPool(SpinupContext_t* ctx, struct Workspace* ws)
{
   if (!ws->isValid)
      destroyWorkspace(ws);
   else {
      ws->isInPool = 1;
      MUTEX_LOCK(ctx->mutex);
      ws->next = ctx->workspaces;
      ctx->workspaces = ws;
      MUTEX_UNLOCK(ctx->mutex);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void cleanVector(PTR vector, struct Workspace* ws)
{
   for (uint32_t j = 0; j < ws->dim; ++j) {
      FEL a = ffExtract(vector, ws->piv[j]);
      if (a != FF_ZERO) {
         PTR basisVector = ffGetPtr(ws->rows, j, ws->noc);
         FEL b = ffExtract(basisVector, ws->piv[j]);
         ffAddMulRow(vector, basisVector, ffNeg(ffDiv(a, b)), ws->noc);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Tries to extend the basis with a "candidate" vector.
/// Returns 0 on success (ws->dim was incremented, rows and pivot table were extended)
/// or -1 if the candidate was a linear combination of the existing basis vectors.
///
/// The caller may store the candidate at the next free workspace position (row index ws->dim + 1)
/// and pass this address in @p vec. This avoids copying the vector internally.
/// 
/// If a spinup script is recorded, «opvec» and «opgen» define the entry to be appended to the
/// script if the cleaned vector is not zero. If no spinup script is recorded «opvec» and «opgen»
/// are ignored.
/// 
/// addToBasis() does not extend the standard basis, this must be done by the caller.

static int addToBasis(struct Workspace* ws, PTR vec, uint32_t opvec, uint32_t opgen)
{
   if (ws->dim >= ws->noc)
      return -1;

   // Copy to next free position if necessary.
   PTR newVec = ffGetPtr(ws->rows, ws->dim, ws->noc);
   if (newVec != vec)
       memcpy(newVec, vec, ffRowSize(ws->noc));

   // Clean with existing basis and add cleaned vector to basis (if not zero).
   cleanVector(newVec, ws);
   if ((ws->piv[ws->dim] = ffFindPivot(newVec, NULL, ws->noc)) == MTX_NVAL) 
      return -1;
   if (ws->ops != NULL) {
      ws->ops[2 * ws->dim] = (int32_t) opvec;
      ws->ops[2 * ws->dim + 1] = (int32_t) opgen;
   }
   ++ws->dim;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Maps the «srcIndex»-th basis vector by all generators in «rep» and tries to extend the subspace
/// basis with the resulting images.
/// If a spinup script is recorded, «opvec» define the source vector value for all script lines
/// that are added.

static void mapAndAddToBasis(struct Workspace* ws, uint32_t srcIndex, const MatRep_t *rep)
{
   PTR srcVec = ffGetPtr(ws->rows, srcIndex, ws->noc);
   for (uint32_t g = 0; g < rep->NGen && ws->dim < ws->noc; ++g) {
      PTR newVec = ffGetPtr(ws->rows, ws->dim, ws->noc);
      ffMapRow(newVec, srcVec, rep->Gen[g]->data, ws->noc, ws->noc);
      const int wasAdded = addToBasis(ws, newVec, srcIndex, g) == 0;
      if (wasAdded && ws->stdRows != NULL) {
         PTR stdSrcVec = ffGetPtr(ws->stdRows, srcIndex, ws->noc);
         PTR stdNewVec = ffGetPtr(ws->stdRows, ws->dim - 1, ws->noc); // dim was already incremented
         ffMapRow(stdNewVec, stdSrcVec, rep->Gen[g]->data, ws->noc, ws->noc);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Transfers the result (subspace + spinup script) from the workspace to the spin up context.
// The workspace becomes invalid and will be deleted by returnWorkspaceToPool().

static void setResult(
      SpinupContext_t* ctx, uint32_t seedVectorNumber, struct Workspace* ws)
{
   ctx->resultSeed = seedVectorNumber;

   Matrix_t* basis = NULL;
   if (ws->stdRows != NULL) {
      basis = matCreateFromBuffer(ws->stdRows, ws->field, ws->dim, ws->noc);
      ws->stdRows = NULL;
   } else {
      basis = matCreateFromBuffer(ws->rows, ws->field, ws->dim, ws->noc);
      ws->rows = NULL;
      matEchelonize(basis);
   }
   if (ctx->submodule != NULL)
      matFree(ctx->submodule);
   ctx->submodule = basis;

   if (ws->ops != NULL) {
      IntMatrix_t* script = imatCreateFromBuffer(ws->ops, ws->dim, 2);
      ws->ops = NULL;
      if (ctx->script != NULL)
         imatFree(ctx->script);
      ctx->script = script;
   }

   ws->isValid = 0;

   MTX_LOG2("Spinup successful for seed #%"PRIu32, seedVectorNumber);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void destroyAllWorkspaces(SpinupContext_t* ctx)
{
    MUTEX_LOCK(ctx->mutex);
    struct Workspace* wsList = ctx->workspaces;
    ctx->workspaces = NULL;
    MUTEX_UNLOCK(ctx->mutex);
    while (wsList != NULL) {
	struct Workspace* ws = wsList;
	wsList = wsList->next;
        ws->isInPool = 0;
        destroyWorkspace(ws);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static struct SpinupContext* createContext(const Matrix_t* seed, const MatRep_t* rep, int flags)
{
   mrValidate(MTX_HERE, rep);
   matValidate(MTX_HERE, seed);
   if (seed->noc == 0) {
      mtxAbort(MTX_HERE, "Dimension is zero");
   if (rep->NGen > 0 && rep->Gen[0]->field != seed->field) {
      mtxAbort(MTX_HERE,
            "Seed is not compatible with representation (field %"PRIu32" vs. %"PRIu32")",
            seed->field, rep->Gen[0]->field);
   }
   }
   if (rep->NGen > 0 && rep->Gen[0]->noc != seed->noc) {
      mtxAbort(MTX_HERE,
            "Seed is not compatible with representation (noc %"PRIu32" vs. %"PRIu32")",
            seed->noc, rep->Gen[0]->noc);
   }

   struct SpinupContext* ctx = ALLOC(SpinupContext_t);
   MUTEX_INIT(ctx->mutex);
   ctx->seed = seed;
   ctx->rep = rep;
   ctx->field = seed->field;
   ctx->noc = seed->noc;
   ctx->flags = flags;
   return ctx;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void destroyContext(struct SpinupContext* ctx)
{
   if (ctx->submodule != NULL) {
      matFree(ctx->submodule);
      ctx->submodule = NULL;
   }
   if (ctx->script != NULL) {
      imatFree(ctx->script);
      ctx->script = NULL;
   }
   destroyAllWorkspaces(ctx);
   MUTEX_DESTROY(ctx->mutex);
   sysFree(ctx);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void spinupOneVector(void* userData, size_t begin, size_t end)
{
   struct SpinupContext* const ctx = (struct SpinupContext*) userData;
   MTX_ASSERT((ctx->flags & SF_MODE_MASK) != SF_COMBINE);
   const uint32_t seedVectorNumber = (uint32_t) begin;
   MTX_ASSERT(seedVectorNumber > 0);

   // If there is already a result for an earlier seed vector, there is nothing to do.
   MUTEX_LOCK(ctx->mutex);
   const int haveOlderResult = (ctx->submodule != NULL && ctx->resultSeed < seedVectorNumber);
   MUTEX_UNLOCK(ctx->mutex);
   if (haveOlderResult)
      return;

   MTX_LOG2("Executing spinup task for seed #%"PRIu32, seedVectorNumber);
   struct Workspace* ws = provideWorkspace(ctx);
   const uint32_t maxSubspaceDim =
      (ctx->maxSubspaceDim > 0) ? ctx->maxSubspaceDim : ws->noc - 1;


   // Start with seed vector
   switch (ctx->flags & SF_SEED_MASK) {
      case SF_FIRST:
      case SF_EACH:
         memcpy(ws->rows, matGetPtr(ctx->seed, seedVectorNumber - 1), ffRowSize(ws->noc));
         break;
      case SF_MAKE:
         svgMake(ws->rows, seedVectorNumber, ctx->seed);
         break;
      default:
         mtxAbort(MTX_HERE, "Unexpected seed mode (flags=0x%04x", (unsigned) ctx->flags);
   }

   // Initialize spin-up with seed vector
   if (ws->stdRows != NULL)
      memcpy(ws->stdRows, ws->rows, ffRowSize(ws->noc));
   if (addToBasis(ws, ws->rows, seedVectorNumber, MTX_NVAL) != 0) {
      mtxAbort(MTX_HERE, "Seed vector is zero");
   }

   // Spin-up loop
   uint32_t srcIndex = 0;
   while (srcIndex < ws->dim && ws->dim <= maxSubspaceDim) {
      mapAndAddToBasis(ws, srcIndex, ctx->rep);
      ++srcIndex;
   }

   int haveResult = 0;
   switch (ctx->flags & SF_MODE_MASK) {
      case SF_SUB:
         haveResult = (srcIndex >= ws->dim && ws->dim <= maxSubspaceDim);
         break;
      case SF_CYCLIC:
         haveResult = (ws->dim == ws->noc);
         break;
      default:
         mtxAbort(MTX_HERE, "Unexpected search mode");
   }
   if (haveResult) {
      // Store the result unless there is already a result for an earlier seed vector.
      MUTEX_LOCK(ctx->mutex);
      if (ctx->submodule == NULL || ctx->resultSeed > seedVectorNumber)
          setResult(ctx, seedVectorNumber, ws);
      MUTEX_UNLOCK(ctx->mutex);
   }
   returnWorkspaceToPool(ctx, ws);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Implementation of seed vector search (spinupFindXxx). Uses threads to spin up multiple seed
// vectors in parallel.
   
static void parallelSpinup(SpinupContext_t* ctx)
{
   ctx->pexGroup = pexCreateGroup();
   
   uint32_t seedVectorNumber = 0;
   int mayAddTasks = 1;
   while (1) {
      // TODO(?): Always wait for the first attempt before producing more tasks?
      // --> mayAddTasks = -1, extend logic in pexThrottle

      // Throttle task creation (keep queue size below 200%)
      pexThrottle(ctx->pexGroup, &mayAddTasks, 200);

      // Get next seed vector number or stop if the seed space is exhausted.
      int haveSeed = 0;
      switch (ctx->flags & SF_SEED_MASK) {
         case SF_FIRST:
            haveSeed = (seedVectorNumber == 0);
            ++seedVectorNumber;
            break;
         case SF_EACH:
            haveSeed = (seedVectorNumber < ctx->seed->nor);
            ++seedVectorNumber;
            break;
         case SF_MAKE:
            haveSeed = svgMakeNext(NULL, &seedVectorNumber, ctx->seed) == 0;
            break;
      }
      if (!haveSeed)
         break;

      // Stop if we have a result
      MUTEX_LOCK(ctx->mutex);
      int hasResult = ctx->submodule != NULL;
      MUTEX_UNLOCK(ctx->mutex);
      if (hasResult) break;

      pexExecuteRange(ctx->pexGroup, spinupOneVector, ctx, seedVectorNumber, 0);
   }

   // Wait for running tasks.
   MTX_LOGD("Stopped at seed #%"PRIu32, seedVectorNumber);
   pexWait(ctx->pexGroup);
   ctx->pexGroup = NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Spin up one or multiple vectors under the action of the generators.

static void simpleSpinup(SpinupContext_t* ctx)
{
   struct Workspace* ws = provideWorkspace(ctx);

   // Add seed vector(s).
   uint32_t nSeedVectors = 0;
   switch (ctx->flags & SF_SEED_MASK) {
      case SF_FIRST: nSeedVectors = 1; break;
      case SF_EACH: nSeedVectors = ctx->seed->nor; break;
      default: mtxAbort(MTX_HERE, "Unsuported seed mode: 0x%04x", ctx->flags);
   }

   for (uint32_t i = 0; i < nSeedVectors; ++i) {
      // Next seed vector
      PTR seedVec = matGetPtr(ctx->seed, i);
      if (ws->stdRows) {
         memcpy(ffGetPtr(ws->stdRows, ws->dim, ws->noc), seedVec, ffRowSize(ws->noc));
      }
      addToBasis(ws, seedVec, i + 1, MTX_NVAL);
      if (ws->dim >= ws->noc)
         break;

      // Spin up
      for (uint32_t srcIndex = 0; srcIndex < ws->dim; ++srcIndex) {
         mapAndAddToBasis(ws, srcIndex, ctx->rep);
      }
      if (ws->dim >= ws->noc)
         break;
   }

   setResult(ctx, 0, ws);
   returnWorkspaceToPool(ctx, ws);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Finds the smallest submodule U of a module V ({0}≤U≤V) containing a given set of vectors.
/// Returns a basis for the submodule.
///
/// @param seed Matrix containing the seed vectors (as rows). The vectors need not be linearly
///    independent.
/// @param rep The matrix representation of the module.

Matrix_t* spinup(const Matrix_t* seed, const MatRep_t* rep)
{
   SpinupContext_t* ctx = createContext(seed, rep, SF_EACH | SF_COMBINE);
   simpleSpinup(ctx);
   Matrix_t* submodule = NULL;
   matMove(&submodule, &ctx->submodule);
   destroyContext(ctx);
   return submodule;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Computes the standard basis for a given set of seed vectors.
///
/// The set S of seed vectors is determined by the rows of @p seed and the value of @p seedMode:
/// - <tt>seedMode = SF_FIRST</tt>: S contains only the first row of @p seed.
/// - <tt>seedMode = SF_EACH</tt>: S contains all rows of @p seed.
/// The rows of @c seed need not be linearly independent.
///
/// The returned matrix is a basis (in standard form, see @ref stdbasis) of the closure of S under
/// the generators. The function does not require that S generates the whole module. If necessary,
/// the caller must check that the returned basis is complete.
///
/// @param seed The seed vector(s).
/// @param rep The matrix representation.
/// @param seedMode Seed mode (SF_FIRST or SF_EACH), see description.
/// @param script Pointer to a variable receiving the spin-up script. May be NULL if the script is
///    not needed.

Matrix_t* spinupStandardBasis(
      IntMatrix_t** script,
      const Matrix_t* seed, const MatRep_t* rep, unsigned seedMode)
{
   if (seed == NULL || seed->nor == 0) {
      mtxAbort(MTX_HERE, "Empty seed");
   }
   switch (seedMode) {
      case SF_FIRST:
      case SF_EACH:
         break;
      default:
         mtxAbort(MTX_HERE, "Invalid seed mode 0x%04x", (unsigned) seedMode);
   }
   SpinupContext_t* ctx = createContext(seed, rep, seedMode | SF_STD);
   simpleSpinup(ctx);
   Matrix_t* basis = NULL;
   matMove(&basis, &ctx->submodule);
   imatMove(script, &ctx->script);
   destroyContext(ctx);
   return basis;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static uint32_t findVector(
      const struct MtxSourceLocation* where,
      Matrix_t** basis,
      const Matrix_t* seed, const MatRep_t* rep, unsigned seedMode, uint32_t maxDimension)
{
   switch (seedMode & SF_SEED_MASK) {
      case SF_FIRST: break;
      case SF_EACH: break;
      case SF_MAKE: break;
      default: mtxAbort(where, "Invalid seed mode 0x%04x", seedMode);
   }
   if (basis == NULL && *basis != NULL) {
      mtxAbort(where, "Basis was not NULL on call");
   }

   SpinupContext_t* ctx = createContext(seed, rep, seedMode);
   if (maxDimension >= ctx->noc)
      mtxAbort(where, "Illegal dimension limit %"PRIu32" (max %"PRIu32")", maxDimension, ctx->noc);
   ctx->maxSubspaceDim = maxDimension;
   parallelSpinup(ctx);
   const uint32_t seedVectorNumber = ctx->resultSeed;
   matMove(basis, &ctx->submodule);
   destroyContext(ctx);
   return seedVectorNumber;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Searches for a seed vector which generates a proper submodule.
/// Returns the seed vector number or 0 if no seed vector generated a proper submodule.
///
/// @param basis
///    If NULL, this parameter is ignored. Otherwise, @p basis must point to a variable which
///    receives the basis created by spinning up the seed vector. The variable must be NULL when
///    spinupFindCyclicVector() is called. If no cyclic vector was found, the variable remains
///    unchanged. The first vector of the returned basis is the seed vector.
/// @param seed The seed vector(s).
/// @param rep The representation.
/// @param seedMode controls how seed vectors are derived from @p seed. The following values
///    can be used:
///    - SF_FIRST
///    - SF_EACH
///    - SF_MAKE
/// @param maxDimension Submodule dimension limit. If not 0, the search is limited to submodules of
///    dimension not greater than @p maxDim. The dimension limit must be strictly less than the
///    module dimension, otherwise the function fails.

uint32_t spinupFindSubmodule(
   Matrix_t** basis,
   const Matrix_t* seed, const MatRep_t* rep, unsigned seedMode, uint32_t maxDimension)
{
   return findVector(MTX_HERE, basis, seed, rep, seedMode | SF_SUB, maxDimension);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Searches for a seed vector which generates the whole module.
/// Returns the seed vector number or 0 if no seed vector generates the whole module.
/// @param basis
///    If NULL, this parameter is ignored. Otherwise, @p basis must point to a variable which
///    receives the basis created by spinning up the seed vector. The variable must be NULL when
///    spinupFindCyclicVector() is called. If no cyclic vector was found, the variable remains
///    unchanged. The first vector of the returned basis is the seed vector.
/// @param seed
///    Seed vector list.
/// @param rep
///    The matrix representation to spin up with.
/// @param seedMode
///    @c SF_FIRST, @c SF_EACH, or @c SF_MAKE. See @ref spinupFindSubmodule.

uint32_t spinupFindCyclicVector(
      Matrix_t** basis,
      const Matrix_t* seed, const MatRep_t* rep, unsigned seedMode)
{
   return findVector(MTX_HERE, basis, seed, rep, seedMode | SF_CYCLIC , 0);
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
