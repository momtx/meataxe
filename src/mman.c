////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Memory Management
////////////////////////////////////////////////////////////////////////////////////////////////////

/// @defgroup mm Memory Management
///
/// This module provides high-level memory management functions such as detecting memory leaks.
///
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <stdlib.h>

/// @private
struct Object {
   struct Object* next;
   struct Object** prev;
   uint32_t seq;
   uint32_t typeId;
};

#if defined(MTX_DEFAULT_THREADS)
static pthread_mutex_t listLock = PTHREAD_MUTEX_INITIALIZER;
#endif
uint32_t sequenceCounter = 0;
static struct Object *objsHead = NULL;
static struct Object **objsTail = &objsHead;
static struct Object *delHead = NULL;
static struct Object **delTail = &delHead;
static size_t nObjs = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Allocates memory for a managed object.
/// The returned memory block has size @a size and starts with an object header (see, for example,
/// @ref Matrix_t) consisting of two pointers and a typeId field. These fields must not be modified
/// by the user.
/// The object must be released with @ref mmFree.

void* mmAlloc(uint32_t typeId, size_t size)
{
    MTX_ASSERT(size >= sizeof(struct Object));
    struct Object* obj = (struct Object*) sysMalloc(size);
    obj->typeId = typeId;
#if defined(MTX_DEFAULT_THREADS)
    pthread_mutex_lock(&listLock);
#endif
    obj->seq = ++sequenceCounter;
    obj->prev = objsTail;
    *objsTail = obj;
    objsTail = &obj->next;
    ++nObjs;
#if defined(MTX_DEFAULT_THREADS)
    pthread_mutex_unlock(&listLock);
#endif
    MTX_LOG2("alloc t=0x%8lx obj=0x%p seq=%lu oc=%lu",
	    (unsigned long) obj->typeId, obj, (unsigned long) obj->seq, (unsigned long) nObjs);
    return obj;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void unlockedMmFree(struct Object* obj, uint32_t typeId)
{
    MTX_ASSERT(nObjs > 0);
    MTX_ASSERT(obj->prev != NULL);
    MTX_ASSERT(obj->typeId == typeId);

    if ((*obj->prev = obj->next) == NULL) {
       if (objsTail == &obj->next)
          objsTail = obj->prev;
       else if (delTail == &obj->next)
          delTail = obj->prev;
       else
          abort();
    } else {
	obj->next->prev = obj->prev;
    }
    --nObjs;

// cannot use logging here!
//    MTX_LOG2("free  t=0x%8lx obj=0x%p seq=%lu oc=%lu", (unsigned long) obj->typeId, obj,
//	    (unsigned long) obj->seq, (unsigned long) nObjs);
    obj->next = NULL;
    obj->prev = NULL;
    obj->typeId = 0;
    sysFree(obj);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Releases object memory.
/// The object passed as first argument must have been created with @ref mmAlloc, using the same
/// type ID.

void mmFree(void* obj, uint32_t typeId)
{
   struct Object*o = (struct Object*) obj;
    #if defined(MTX_DEFAULT_THREADS)
    pthread_mutex_lock(&listLock);
    #endif
    unlockedMmFree(o, typeId);
    #if defined(MTX_DEFAULT_THREADS)
    pthread_mutex_unlock(&listLock);
    #endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Prints a warning is there are any live (allocated but not released) objects.

void mmLeakCheck()
{
   if (objsHead == NULL) {
      return;
   }
   size_t nMatrix = 0;
   size_t nPolynomial = 0;
   size_t nFPoly = 0;
   size_t nPermutation = 0;
   size_t nIntMatrix = 0;
   size_t nBinFile = 0;
   size_t nBsFixed = 0;
   size_t nBsDynamic = 0;
   size_t nWgen = 0;
   size_t nStFile = 0;
   size_t nOther = 0;
   size_t nCharpol = 0;

   for (const struct Object* obj = objsHead; obj != NULL; obj = obj->next) {
      switch (obj->typeId) {
         case MTX_TYPE_MATRIX: ++nMatrix; break;
         case MTX_TYPE_PERMUTATION: ++nPermutation; break;
         case MTX_TYPE_POLYNOMIAL: ++nPolynomial; break;
         case MTX_TYPE_FPOLY: ++nFPoly; break;
         case MTX_TYPE_INTMATRIX: ++nIntMatrix; break;
         case MTX_TYPE_BINFILE: ++nBinFile; break;
         case MTX_TYPE_STFILE: ++nStFile; break;
         case MTX_TYPE_WORD_GENERATOR: ++nWgen; break;
         case MTX_TYPE_BITSTRING_FIXED: ++nBsFixed; break;
         case MTX_TYPE_BITSTRING_DYNAMIC: ++nBsDynamic; break;
         case MTX_TYPE_CPSTATE: ++nCharpol; break;
         default: ++nOther; break;
      }
      MTX_LOGE("leak t=0x%8lx obj=0x%p seq=%lu",
	    (unsigned long) obj->typeId, obj, (unsigned long) obj->seq);
   }
   MTX_XLOGE(msg) {
      sbAppend(msg, "Leak check:");
      if (nMatrix > 0) sbPrintf(msg," %lu mat", (unsigned long)nMatrix);
      if (nPermutation > 0) sbPrintf(msg," %lu perm", (unsigned long)nPermutation);
      if (nPolynomial > 0) sbPrintf(msg," %lu pol", (unsigned long)nPolynomial);
      if (nFPoly > 0) sbPrintf(msg," %lu fpol", (unsigned long)nFPoly);
      if (nIntMatrix > 0) sbPrintf(msg," %lu imat", (unsigned long)nIntMatrix);
      if (nBinFile > 0) sbPrintf(msg," %lu binfile", (unsigned long)nBinFile);
      if (nStFile > 0) sbPrintf(msg," %lu stfile", (unsigned long)nBinFile);
      if (nWgen > 0) sbPrintf(msg," %lu wg", (unsigned long)nWgen);
      if (nBsFixed > 0) sbPrintf(msg," %lu bsfix", (unsigned long)nBsFixed);
      if (nBsDynamic > 0) sbPrintf(msg," %lu bsdyn", (unsigned long)nBsDynamic);
      if (nCharpol > 0) sbPrintf(msg," %lu chpol", (unsigned long)nCharpol);
      if (nOther > 0) sbPrintf(msg," %lu other", (unsigned long)nOther);
   }
   mtxAbort(MTX_HERE, "leak check failed");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

//static struct Object* mmTail()
//{
//   if (objsTail == &objsHead) return NULL;
//   return (struct Object*) objsTail; // this works because "next" is the first field in struct Object!
//}

////////////////////////////////////////////////////////////////////////////////////////////////////

uint32_t mmCheckpoint()
{
#if defined(MTX_DEFAULT_THREADS)
    pthread_mutex_lock(&listLock);
#endif
    uint32_t checkpoint = sequenceCounter;
#if defined(MTX_DEFAULT_THREADS)
    pthread_mutex_unlock(&listLock);
#endif
   return checkpoint;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void destroy(void* object)
{
   struct Object* obj = (struct Object*)object;
   switch (obj->typeId) {
      case MTX_TYPE_BINFILE: mfClose((MtxFile_t*)obj);
         return;
      case MTX_TYPE_BITSTRING_FIXED:
      case MTX_TYPE_BITSTRING_DYNAMIC: bsFree((BitString_t*)obj); return;
      case MTX_TYPE_CPSTATE: charpolFree((Charpol_t*) obj); return;
      case MTX_TYPE_INTMATRIX: imatFree((IntMatrix_t*) obj); return;
      case MTX_TYPE_MATREP: mrFree((MatRep_t*) obj); return;
      case MTX_TYPE_MATRIX: matFree((Matrix_t*) obj); return;
      case MTX_TYPE_PERMUTATION: permFree((Perm_t*) obj); return;
      case MTX_TYPE_POLYNOMIAL: polFree((Poly_t*) obj); return;
      case MTX_TYPE_FPOLY: fpFree((FPoly_t*) obj); return;
      case MTX_TYPE_STFILE: stfClose((StfData*) obj); return;
      case MTX_TYPE_STRBUF: sbFree((StrBuffer_t*) obj); return;
      case MTX_TYPE_WORD_GENERATOR: wgFree((WgData_t*)obj); return;
   }
   abort();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Destroys all recent objects up to, but not including, @a checkpoint.
/// If @a checkpoint was already destroyed, nothing happens.

void mmRollback(uint32_t checkpoint)
{
#if defined(MTX_DEFAULT_THREADS)
    pthread_mutex_lock(&listLock);
#endif
    if (delHead != NULL || delTail != &delHead)
       abort();
    struct Object** o = &objsHead;
    while (*o != NULL && (*o)->seq <= checkpoint)
       o = &(*o)->next;
    if (*o != NULL && (*o)->seq > checkpoint) {
       // move to delete list
       delHead = *o;
       (*o)->prev = &delHead;
       delTail = objsTail;
       objsTail = o;
       *objsTail = NULL;
    }
#if defined(MTX_DEFAULT_THREADS)
    pthread_mutex_unlock(&listLock);
#endif

    while (delHead != NULL) {
       destroy(delHead);
    }
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
