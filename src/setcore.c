////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Basic set functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


static const int InitialSize = 10;
static const unsigned long SetMagic = 0xEF452338;

/// @defgroup intset Sets of Integers.
/// @{

/// @class Set_t
/// @brief
/// A Set of Integers.
/// The Set_t structure represents a set of (long) integers. Internally, the set is stored
/// as a sorted list. Insert operations are relatively expensive, especially for large sets.
/// So, if you expect a lot of inserts, the BitString_t data type may be a better choice.

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns 1 if the given set is valid and 0 otherwise,

int setIsValid(const Set_t *s)
{
   return s != NULL && s->Magic == SetMagic && s->Size >= 0
      && s->BufSize >= s->Size && s->Data != NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Checks whether the given set is valid and aborts the program if the test fails.

void setValidate(const struct MtxSourceLocation* src, const Set_t *s)
{
   if (s == NULL) {
      mtxAbort(MTX_HERE,"NULL set");
   }
   if ((s->Magic != SetMagic) || (s->Size < 0) || (s->BufSize < s->Size)) {
      mtxAbort(MTX_HERE,"Invalid set (Magic=%d, Size=%d, BufSize=%d)",
                 (int)s->Magic,s->Size,s->BufSize);
   }
   if (s->Data == NULL) {
      mtxAbort(MTX_HERE,"Data==NULL in set");
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Create a new set.
/// This function creates a new, empty set. To destroy a set,
/// use setFree(), @em not sysFree().
/// @return Pointer to the new set or 0 on error.

Set_t *setAlloc()
{
   Set_t *x;

   x = ALLOC(Set_t);
   if (x == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate set");
      return NULL;
   }
   x->Size = 0;
   x->BufSize = InitialSize;
   x->Data = NALLOC(long,InitialSize);
   if (x->Data == NULL) {
      sysFree(x);
      mtxAbort(MTX_HERE,"Cannot allocate set data");
      return NULL;
   }
   x->Magic = SetMagic;
   return x;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Destroy a set.
/// This function frees an integer set. The argument must be a Set_t
/// structure which has previously been allocated with setAlloc().
/// @param x Pointer to the set.
/// @return 0 on success, -1 on error.

int setFree(Set_t *x)
{
   setValidate(MTX_HERE, x);
   sysFree(x->Data);
   memset(x,0,sizeof(Set_t));
   sysFree(x);
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Duplicate a set.
/// @param s Pointer to the set.
/// @return Pointer to a copy of the set, or 0 on error.

Set_t *setDup(const Set_t *s)
{
   Set_t *x;

   setValidate(MTX_HERE, s);
   x = ALLOC(Set_t);
   if (x == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate set");
      return NULL;
   }

   x->Size = s->Size;
   x->BufSize = s->Size;
   x->Data = NALLOC(long,x->BufSize);
   if (x->Data == NULL) {
      sysFree(x);
      mtxAbort(MTX_HERE,"Cannot allocate set data");
      return NULL;
   }
   memcpy(x->Data,s->Data,sizeof(x->Data[0]) * s->Size);
   x->Magic = SetMagic;
   return x;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
