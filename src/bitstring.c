////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Bit Strings, core functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

/// @defgroup bs Bit Strings
/// @{
/// @details
/// The @ref BitString_t type represents a string of 0's and 1's. You may also think of it as a set
/// of nonnegative integers, where each 1 in the bit string means that the set contains the
/// corresponding number.
///
/// Variable bit strings are with @ref bsAllocEmpty. They are initially empty and grow as needed
/// when bits are set. The @ref bsTrim function can be used to free unused memory occupied by
/// trailing 0 bits.
///
/// Alternativley, a fixed-size bit string can be created using @ref bsAlloc. This is useful if the
/// number of bits is known, e.g., when dealing with subsets of a given finite set M. Fixed-size
/// bit strings do not allow access to bits beyond the fixed range 0...N-1 and will abort the
/// program instead of extending the bit string.

////////////////////////////////////////////////////////////////////////////////////////////////////

static const size_t BPL = sizeof(long) * 8;

int bsIsValid(const BitString_t *bs)
{
   return bs != NULL
      && (bs->typeId == MTX_TYPE_BITSTRING_FIXED || bs->typeId == MTX_TYPE_BITSTRING_DYNAMIC)
      && bs->capacity >= bs->size;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Checks if a bit string is valid and aborts the program if the test fails.

void bsValidate(const struct MtxSourceLocation* src, const BitString_t *bs)
{
   if (!bsIsValid(bs))
      mtxAbort(src,"invalid bit string");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Creates an empty, variable-sized bit string. The initial bit string length is 0.

BitString_t *bsAllocEmpty()
{
   BitString_t *bs = (BitString_t *) mmAlloc(MTX_TYPE_BITSTRING_DYNAMIC, sizeof(BitString_t));
   bs->size = 0;
   bs->capacity = 0;
   bs->data = NULL;
   return bs;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Creates a fixed-size bit string.
/// All bits are initially cleared. The bit string size cannot be changed. Accessing bits at
/// position @a size or greater will result in an error and abort the program.
///
/// See also @ref bsFree

BitString_t *bsAlloc(size_t size)
{
   BitString_t *bs = (BitString_t *) mmAlloc(MTX_TYPE_BITSTRING_FIXED, sizeof(BitString_t));
   bs->size = size;
   bs->capacity = sysPad(size, BPL);
   bs->data = NALLOC(uint8_t, bs->capacity / 8);
   return bs;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Destroys a bit string and releases all associated memory.

int bsFree(BitString_t *bs)
{
   bsValidate(MTX_HERE, bs);
   sysFree(bs->data);
   bs->data = NULL;
   mmFree(bs, bs->typeId);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void badIndex(const struct MtxSourceLocation* src, const BitString_t* bs, size_t i)
{
   mtxAbort(MTX_HERE, "Bit string index out of range: i=%lu size=%lu",
         (unsigned long)i, (unsigned long)bs->size);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the value of a bit in a bit string.

int bsTest(const BitString_t *bs, size_t i)
{
   bsValidate(MTX_HERE, bs);
   if (bs->typeId == MTX_TYPE_BITSTRING_FIXED) {
      if (i >= bs->size)
         badIndex(MTX_HERE, bs, i);
   } else {
      if (i >= bs->capacity)
         return 0;
   }
   return (bs->data[i / 8] & (0x80 >> (i % 8))) != 0 ? 1 : 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void resizeBuffer(BitString_t* bs, size_t newCapacity)
{
   MTX_ASSERT(newCapacity % BPL == 0);
   bs->data = NREALLOC(bs->data, uint8_t, newCapacity / 8);
   if (newCapacity > bs->capacity) {
      memset(bs->data + bs->capacity/8, 0, (newCapacity - bs->capacity) / 8);
   }
   bs->capacity = newCapacity;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Free unused memory in a bit string.
/// This function has no effect on fixed-size bit strings.
/// It may be used on variable-sized bit strings after clearing bits to release memory occupied by
/// trailing 0 bits.

void bsTrim(const BitString_t* bs)
{
   bsValidate(MTX_HERE, bs);
   if (bs->typeId == MTX_TYPE_BITSTRING_FIXED)
      return;
   unsigned long* lData = (unsigned long*) bs->data;
   size_t lSize = bs->capacity / BPL;
   while (lSize > 0 && lData[lSize - 1] == 0)
      --lSize;
   size_t newCapacity = lSize * BPL;
   if (newCapacity < bs->capacity) {
      resizeBuffer((BitString_t*)bs, newCapacity);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Changes and locks the size of a bit string.
/// This function can be called on any bit string, variable or fixed size.
/// It converts the bit string into a fixed-size bit string of the given size. Existing bits at
/// position @a newSize and above will be lost. If the size is increased, new bits are initialized
/// with 0.

void bsResize(BitString_t* bs, size_t newSize)
{
   bsValidate(MTX_HERE, bs);

   // Extend/shrink the buffer.
   const size_t newCapacity = sysPad(newSize, BPL);
   resizeBuffer(bs, newCapacity);

   // Clear insignificant bits.
   for (size_t i = newSize + 1; i < newCapacity; ++i) {
      bs->data[i / 8] &= ~(0x80 >> (i % 8));
   }
   bs->size = newSize;
   bs->typeId = MTX_TYPE_BITSTRING_FIXED;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Sets a bit in a bit string.
/// 
/// If the bit string is of fixed size, @a i must be less than <tt>bs->size</tt>.

void bsSet(BitString_t* bs, size_t i)
{
   bsValidate(MTX_HERE, bs);
   if (bs->typeId == MTX_TYPE_BITSTRING_FIXED && i >= bs->size) {
      badIndex(MTX_HERE, bs, i);
   }
   if (i >= bs->capacity) {
      resizeBuffer(bs, sysPad(i + 1, BPL));
   }
   bs->data[i / 8] |= 0x80 >> (i % 8);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Clears a bit in a bit string.
/// 
/// If the bit string is of fixed size, @a i must be less than <tt>bs->size</tt>.

void bsClear(BitString_t* bs, size_t i)
{
   bsValidate(MTX_HERE, bs);
   if (bs->typeId == MTX_TYPE_BITSTRING_FIXED && i >= bs->size) {
      badIndex(MTX_HERE, bs, i);
   }
   if (i < bs->capacity)
      bs->data[i / 8] &= ~(0x80 >> (i % 8));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Clears all bits in a bit string.
/// If the bit string is dynamic, its size is reset to zero.
/// @return 0 on success, -1 on error.

void bsClearAll(BitString_t *bs)
{
   bsValidate(MTX_HERE, bs);
   if (bs->typeId == MTX_TYPE_BITSTRING_DYNAMIC) {
      bs->size = 0;
      bs->capacity = 0;
      sysFree(bs->data);
      bs->data = NULL;
   } else {
      memset(bs->data,0,bs->capacity / 8);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void validate2(const BitString_t* a, const BitString_t* b)
{
   bsValidate(MTX_HERE, a);
   bsValidate(MTX_HERE, b);
   if (a->typeId != b->typeId)
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
   if (a->typeId == MTX_TYPE_BITSTRING_FIXED && a->size != b->size)
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Performs a logical and (intersection) on two bit strings and stores the result in @a dest.
/// @a dest and @a src must either both be of fixed-size type or both be dynamic. If they are of
/// fixed-size type, their size must be equal. The type of @a dest remains unchanged.

void bsAnd(BitString_t* dest, const BitString_t* src)
{
   validate2(dest, src);
   if (dest->typeId == MTX_TYPE_BITSTRING_DYNAMIC) {
      // variable type - shrink dest if possible
      if (dest->capacity > src->capacity) {
         resizeBuffer(dest, src->capacity);
      }
   }
   unsigned long *dp = (unsigned long*)dest->data;
   const unsigned long *sp = (const unsigned long*)src->data;
   for (size_t i = dest->capacity / BPL; i > 0; --i)
      *dp++ &= *sp++;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Performs a logical or (union) on two bit strings and stores the result in @a dest.
/// @a dest and @a src must either both be of fixed-size type or both be dynamic. If they are of
/// fixed-size type, their size must be equal. The type of @a dest remains unchanged.

void bsOr(BitString_t* dest, const BitString_t* src)
{
   validate2(dest, src);

   if (dest->typeId == MTX_TYPE_BITSTRING_DYNAMIC) {
      // variable type - extend dest if necessary
      if (dest->capacity < src->capacity) {
         resizeBuffer(dest, src->capacity);
      }
   }
   unsigned long *dp = (unsigned long*)dest->data;
   const unsigned long *sp = (const unsigned long*)src->data;
   for (size_t i = dest->capacity / BPL; i > 0; --i)
      *dp++ |= *sp++;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Clears all bits in @a dest which are also set in @a src.
/// @a dest and @a src must either both be of fixed-size type or both be dynamic. If they are of
/// fixed-size type, their size must be equal. The type of @a dest remains unchanged.

void bsMinus(BitString_t *dest, const BitString_t *src)
{
   validate2(dest,src);

   const size_t minCapacity = src->capacity < dest->capacity ? src->capacity : dest->capacity;
   unsigned long *dp = (unsigned long*)dest->data;
   const unsigned long *sp = (const unsigned long*)src->data;
   for (size_t i = minCapacity / BPL; i > 0; --i)
      *dp++ &= ~*sp++;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns 1 if @a a is a subset of @a b and 0 otherwise.
///
/// @a dest and @a src must either both be of fixed-size type or both be dynamic. If they are of
/// fixed-size type, their size must be equal.

int bsIsSub(const BitString_t *a, const BitString_t *b)
{
   validate2(a, b);

   // Check common range.
   const unsigned long *ap = (const unsigned long*)a->data;
   const unsigned long *bp = (const unsigned long*)b->data;
   const size_t minCapacity = a->capacity < b->capacity ? a->capacity : b->capacity;
   for (size_t i = minCapacity / BPL; i > 0; --i) {
      if (*ap != (*ap & *bp))
         return 0;
      ++ap;
      ++bp;
   }

   // For variable bit strings, a may be longer than b.
   if (a->capacity > minCapacity) {
      for (size_t i = (a->capacity - minCapacity) / BPL; i > 0; --i) {
         if (*ap++ != 0)
            return 0;
      }
   }

   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Counts matching bits in two bit strings.
///
/// This function calculates the cardinality of the intersection of two bit strings, i.e.,
/// the number of bits that are set in both @a a and @a b.
/// @a dest and @a src must either both be of fixed-size type or both be dynamic. If they are of
/// fixed-size type, their size must be equal.

size_t bsIntersectionCount(const BitString_t *a, const BitString_t *b)
{
   validate2(a, b);

   static const int BitCount[256] =  /* Number of '1' bits in binary representation */
   {
      0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
      1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
      1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
      2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
      1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
      2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
      2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
      3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
      1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
      2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
      2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
      3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
      2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
      3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
      3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
      4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
   };
   size_t count = 0;

   const unsigned long *ap = (const unsigned long*)a->data;
   const unsigned long *bp = (const unsigned long*)b->data;
   const size_t minCapacity = a->capacity < b->capacity ? a->capacity : b->capacity;
   for (size_t i = minCapacity / BPL; i > 0; --i) {
      for (unsigned long x = *ap++ & *bp++; x != 0; x >>= 8) {
         count += BitCount[x & 0xFF];
      }
   }
   return count;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Compare two bit strings.
/// @a dest and @a src must either both be of fixed-size type or both be dynamic.
/// @return 0 if the bit strings are equal. Otherwise the return value is different from zero.

int bsCompare(const BitString_t *a, const BitString_t *b)
{
   bsValidate(MTX_HERE, a);
   bsValidate(MTX_HERE, b);
   if (a->typeId != b->typeId)
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);

   const unsigned long *ap = (const unsigned long*)a->data;
   const unsigned long *bp = (const unsigned long*)b->data;
   const size_t minCapacity = a->capacity < b->capacity ? a->capacity : b->capacity;
   for (size_t i = minCapacity / BPL; i > 0; --i) {
      if (*ap != *bp) {
         // Compare bytewise to make the result independent of the endianness.
         const uint8_t *a8 = (const uint8_t *)ap;
         const uint8_t *b8 = (const uint8_t *)bp;
         for (unsigned k = sizeof(long); k > 0; --k) {
            if (*a8 < *b8) return -1;
            if (*a8 > *b8) return 1;
            ++a8;
            ++b8;
         }
         MTX_ASSERT(0); // cannot happen
      }
      ++ap;
      ++bp;
   }
   if (a->capacity > minCapacity) {
      for (size_t i = (a->capacity - minCapacity) / BPL; i > 0; --i) {
         if (*ap++ != 0)
            return 1;
      }
   }
   else if (b->capacity > minCapacity) {
      for (size_t i = (b->capacity - minCapacity) / BPL; i > 0; --i) {
         if (*bp++ != 0)
            return -1;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Copies a bit string.
/// The current contents of @a dest is replaced with a copy of @a src. This includes the bit string
/// type (variable or fixed-size).

void bsCopy(BitString_t *dest, const BitString_t *src)
{
   bsValidate(MTX_HERE, src);
   bsValidate(MTX_HERE, dest);
   dest->typeId = src->typeId;
   dest->size = src->size;
   dest->capacity = src->capacity;
   dest->data = NREALLOC(dest->data, uint8_t, dest->capacity / 8);
   memcpy(dest->data, src->data, dest->capacity / 8);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Creates a copy of a bit string.
/// The copy has the same type (fixed-size or dynamic) as the source.

BitString_t *bsDup(const BitString_t *src)
{
   bsValidate(MTX_HERE, src);
   BitString_t *n = NULL;
   if (src->typeId == MTX_TYPE_BITSTRING_FIXED)
      n = bsAlloc(src->size);
   else {
      n = bsAllocEmpty();
      resizeBuffer(n, src->capacity);
   }
   memcpy(n->data,src->data,n->capacity / 8);
   return n;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Prints a bit string on stdout.

void bsPrint(const char *name, const BitString_t *bs)
{
   if (name != NULL) {
      printf("%s=",name);
   }
   size_t end = bs->size;
   if (bs->typeId == MTX_TYPE_BITSTRING_DYNAMIC) {
      uint8_t *b = bs->data + bs->capacity / 8;
      while (b > bs->data && b[-1] == 0) --b;
      end = (b - bs->data) * 8;
   } 

   for (size_t i = 0; i < end; ++i) {
      putc(bsTest(bs,i) ? '1' : '0', stdout);
   }
   putc('\n', stdout);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes a bit string to a file.
/// The file must be open for writing.

void bsWrite(const BitString_t *bs, MtxFile_t* file)
{
   bsValidate(MTX_HERE, bs);
   MTX_ASSERT(file != NULL);

   uint32_t fileHeader[3] = {0,0,0};
   if (bs->typeId == MTX_TYPE_BITSTRING_FIXED) {
      fileHeader[0] = MTX_TYPE_BITSTRING_FIXED;
      fileHeader[1] = bs->size;
   } else {
      bsTrim(bs);
      fileHeader[0] = MTX_TYPE_BITSTRING_DYNAMIC;
      fileHeader[1] = bs->capacity;
   }
   mfWrite32(file, fileHeader, 3);
   mfWrite8(file, bs->data, sysPad(fileHeader[1], 8) / 8);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static uint32_t checkType(MtxFile_t *f)
{
   const uint32_t objectType = mfObjectType(f);
   if (objectType != MTX_TYPE_BITSTRING_FIXED && objectType != MTX_TYPE_BITSTRING_DYNAMIC) {
      mtxAbort(MTX_HERE, "%s: bad type 0x%lx, expected 0x%lx or 0x%lx (BITSTRING)",
         f->name, (unsigned long) objectType,
         (unsigned long) MTX_TYPE_BITSTRING_FIXED,
         (unsigned long) MTX_TYPE_BITSTRING_DYNAMIC);
   }
   return objectType;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

BitString_t *bsReadData(MtxFile_t *f)
{
   const uint32_t objectType = checkType(f);
   BitString_t *bs = bsAlloc(f->header[1]);
   if (objectType == MTX_TYPE_BITSTRING_DYNAMIC) {
      bs->typeId = MTX_TYPE_BITSTRING_DYNAMIC;
   }
   mfRead8(f, bs->data, sysPad(bs->size, 8) / 8);

   // Make sure a second read attempt will fail.
   f->header[0] = 0xFFFFFFFF;

   return bs;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads a bit string from a file and returns the bit string.
/// The file must be open for reading.

BitString_t* bsRead(MtxFile_t* f)
{
   mfReadHeader(f);
   return bsReadData(f);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void bsSkip(MtxFile_t* f)
{
   mfReadHeader(f);
   checkType(f);
   const size_t bsSize = f->header[1];
   const size_t nBytes = sysPad(bsSize, 8) / 8;
   sysFseekRelative(f->file, nBytes);
   f->header[0] = 0xFFFFFFFF;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Finds the first "1" bit in a bit string.
/// If the bit string contains at least one "1" bit, the index of the first "1" is stored in the
/// variable pointed to by @a indexVar and the function returns 1. Otherwise the function returns 0
/// and the variable pointed to by @a indexVar is not changed.

int bsFirst(const BitString_t* bs, size_t* indexVar)
{
   const unsigned long* lp = (const unsigned long*) bs->data;
   const unsigned long* const end = (bs->typeId == MTX_TYPE_BITSTRING_FIXED) ?
      lp + sysPad(bs->size, BPL) / BPL : lp + bs->capacity / BPL;
   while (lp < end && *lp == 0) ++lp;
   if (lp >= end)
      return 0;

   size_t i = (lp - (const unsigned long*) bs->data) * BPL;
   for (const uint8_t* bp = (const uint8_t*) lp; ; ++bp, i += 8) {
      if (*bp & 0x80) { *indexVar = i; break; }
      if (*bp & 0x40) { *indexVar = i + 1; break; }
      if (*bp & 0x20) { *indexVar = i + 2; break; }
      if (*bp & 0x10) { *indexVar = i + 3; break; }
      if (*bp & 0x08) { *indexVar = i + 4; break; }
      if (*bp & 0x04) { *indexVar = i + 5; break; }
      if (*bp & 0x02) { *indexVar = i + 6; break; }
      if (*bp & 0x01) { *indexVar = i + 7; break; }
   }
   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Finds the next "1" bit starting at a given position. 

/// If @em s is the value of the "index variable" pointed to by @a indexVar, the function finds the
/// first "1" at a position greater than @em s, stores its position in the index variable and
/// returns 1. If no such bit is found, the function returns 0 and does not change the index
/// variable.

int bsNext(const BitString_t* bs, size_t* indexVar)
{
   const uint8_t* const end = bs->data +
      ((bs->typeId == MTX_TYPE_BITSTRING_FIXED) ? sysPad(bs->size, 8) / 8 : bs->capacity / 8);

   // Start with the bit following *indexVar.
   size_t i = *indexVar + 1;
   const uint8_t* bp = bs->data +  i / 8;

   // Check if a bit is set in the same byte.
   if (i % 8 != 0) {
      while (i % 8 != 0) {
         if (*bp & (0x80 >> (i % 8))) {
            *indexVar = i;
            return 1;
         }
         ++i;
      }
      ++bp;
   }

   // Search next "1" bit in the remaining data.
   while (bp < end) {
      if (i % BPL == 0 && *(const unsigned long*)bp == 0) {
         bp += sizeof(long);
         i += BPL;
      }
      else if (*bp == 0) {
         ++bp;
         i += 8;
      } else {
         if (*bp & 0x80) { *indexVar = i; return 1; }
         if (*bp & 0x40) { *indexVar = i + 1; return 1; }
         if (*bp & 0x20) { *indexVar = i + 2; return 1; }
         if (*bp & 0x10) { *indexVar = i + 3; return 1; }
         if (*bp & 0x08) { *indexVar = i + 4; return 1; }
         if (*bp & 0x04) { *indexVar = i + 5; return 1; }
         if (*bp & 0x02) { *indexVar = i + 6; return 1; }
         if (*bp & 0x01) { *indexVar = i + 7; return 1; }
         MTX_ASSERT(0);
      }
   }

   // Nothing found
   return 0;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
