////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Permutations
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

/// @defgroup perm Permutations
/// @details
/// In the MeatAxe, a permutation of degree n operates on {0,1,...,n-1} and is represented
/// by a Perm_t structure.
/// However, in the textual representation produced by permPrint() or by the
/// @ref prog_zpr "zpr" program, the points are numbered from 1...n.
///
/// Only permutations of equal degree can be multiplied. This can be confusing
/// because the textual representation produced by permPrint() does not include the
/// degree, and fixed points are always suppressed. For example "(3 4)(5 6 7)" could
/// be interpreted as a permutation of degree 8 or any higher degree. All these permutations
/// are in some natural way equal to each other, but they are different and incompatible
/// in the MeatAxe.
///
/// Permutations are usually created with permAlloc() or read from a file with permRead().
/// When a permutation is no longer used, the application must release the associated memory
/// by calling permFree().
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// @class Perm_t
/// @details
/// Internally, a permutation is represented as an array of 32-bit integers containing the
/// images of 0,1,...,n-1. The maximum degree is 2^32-1.

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Compare two permutations.
///
/// This function returns -1, 0, or 1 if the permutation @a a is less that, equal to, or greater
/// than permutation @a b.  The ordering of permutations is defined as follows:
/// - If the permutations have different degrees, the permutations with the smaller degree is
///   smaller.
/// - Otherwise, the order is defined by the lexicograhical order of the sequences
///   (a(0), … , a(n-1)) and (b(0),…,b(n-1)), where n is the degree.

int permCompare(const Perm_t *a, const Perm_t *b)
{
   // check arguments
   permValidate(MTX_HERE,a);
   permValidate(MTX_HERE,b);

   // compare degrees
   if (a->degree > b->degree) { return 1; }
   if (a->degree < b->degree) { return -1; }
  
   // compare the entries
   uint32_t* pa = a->data;
   uint32_t* pb = b->data;
   while (pa < a->data + a->degree) {
      if (*pa > *pb) return 1;
      if (*pa < *pb) return -1;
      ++pa;
      ++pb;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check a permutation.
/// This function checks if the argument is a pointer to a valid permutation.
/// If the permutation is o.k., the function returns 1.
/// Otherwise, an error is signalled and, if the error handler does not
/// terminate the program, the function returns 0.
/// @return 1 if @a p is a valid permutation, 0 otherwise.

int permIsValid(const Perm_t *p)
{
   if (p == NULL) {
      return 0;
   }
   if (p->typeId != MTX_TYPE_PERMUTATION || p->degree < 0 || p->data == NULL) {
      return 0;
   }
   for (int i = 0; i < p->degree; ++i) {
      if (p->data[i] < 0 || p->data[i] >= p->degree) {
         return 0;
      }
   }

   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Aborts the program if the passed permutation object is not valid.

void permValidate(const struct MtxSourceLocation* src, const Perm_t* p)
{
   if (p == NULL) {
      mtxAbort(src, "NULL permutation");
   }
   if (p->typeId != MTX_TYPE_PERMUTATION || p->data == NULL) {
      mtxAbort(src, "Invalid permutation (type=0x%lu, deg=%lu)",
         (unsigned long) p->typeId,
         (unsigned long) p->degree);
   }
   for (int i = 0; i < p->degree; ++i) {
      if (p->data[i] < 0 || p->data[i] >= p->degree) {
         mtxAbort(src, "Invalid value %d in permutation (deg = %d)", (int) p->data[i], p->degree);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Allocate a permutation
/// This function creates a permutation of the specified degree.
/// The new permutation is initialized to the identity.
/// @param deg Degree.
/// @return Pointer to the new permutation, or 0 on error.

Perm_t *permAlloc(uint32_t deg)
{
   Perm_t* p = (Perm_t*) mmAlloc(MTX_TYPE_PERMUTATION, sizeof(Perm_t));
   p->degree = deg;
   p->data = NALLOC(uint32_t,deg);
   for (uint32_t i = 0; i < deg; ++i) {
      p->data[i] = i;
   }
   return p;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Free a permutation.
/// This function deletes a permutation and returns the memory to the
/// system. Do not use sysFree() on permutations because this would only
/// free the Perm_t structure but not the data buffer.
/// @param p Pointer to the permutation.
/// @return 0 on success, -1 on error.

void permFree(Perm_t *p)
{
   permValidate(MTX_HERE, p);

   sysFree(p->data);
   p->data = NULL;
   mmFree(p, MTX_TYPE_PERMUTATION);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Duplicate a permutation.
/// This function creates a copy of an existing permutation.
/// @param src Pointer to the permutation.
/// @return Pointer to a copy of @em src or 0 on error.

Perm_t *permDup(const Perm_t *src)
{
   permValidate(MTX_HERE, src);
   Perm_t *p = permAlloc(src->degree);
   if (p == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate result");
      return NULL;
   }
   memcpy(p->data,src->data,(size_t) src->degree * sizeof(p->data[0]));
   return p;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Inverse of a permutation
/// This function calulates the inverse of a permutation.
/// @param src Pointer to the permutation.
/// @return The inverse of @em src, or 0 on error.

Perm_t *permInverse(const Perm_t *src)
{
   permValidate(MTX_HERE, src);

   Perm_t* inv = permAlloc(src->degree);
   uint32_t* d = inv->data;
   const uint32_t* s = src->data + src->degree;
   for (uint32_t i = src->degree; i > 0; ) {
      --i;
      --s;
      d[*s] = i;
   }

   return inv;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Multiply permutations.
/// This function multiplies @em dest from the right by @em src. Both
/// permutations must have the same degree.
/// @param dest Pointer to the first permutation.
/// @param src Pointer to the second permutation.
/// @return @em dest, or 0 on error.

Perm_t *permMul(Perm_t *dest, const Perm_t *src)
{
   permValidate(MTX_HERE, src);
   permValidate(MTX_HERE, dest);
   if (dest->degree != src->degree) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_INCOMPAT);
   }

   uint32_t* d = dest->data;
   const uint32_t* const s = src->data;
   for (uint32_t i = dest->degree; i > 0; --i) {
      *d = s[*d];
      ++d;
   }
   return dest;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the order of a permutation.
///
/// The function may fail if the order is too big.

uint32_t permOrder(const Perm_t *perm)
{

   permValidate(MTX_HERE, perm);
   if (perm->degree < 2) {
      return 1;
   }

   const uint32_t deg = perm->degree;
   const uint32_t* p = perm->data;

   uint8_t* done = NALLOC(uint8_t, deg);
   uint8_t* seed = done - 1;
   uint8_t* seedEnd = done + deg;

   // Calculate the order by running through all orbits. 
   uint32_t order = 1;
   while (1) {
      // Find start of next orbit.
      while (*++seed != 0 && seed != seedEnd);
      if (seed == seedEnd)
         break;             // Done!
      
      // Find orbit size
      uint32_t orbitSize = 0;
      uint32_t x = (uint32_t)(seed - done);
      while (1) {
         MTX_ASSERT_DEBUG(x >= 0 && x < deg);
         if (done[x])
            break;      // orbit complete
         done[x] = 1;
         x = p[x];
         ++orbitSize;
      }

      // Calculate the l.c.m of all orbit sizes
      order = lcm32u(order,orbitSize);
   }

   sysFree(done);

   return order;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int width(uint32_t x)
{
   if (x < 10) return 1;
   if (x < 100) return 2;
   if (x < 1000) return 3;
   if (x < 10000) return 4;
   if (x < 100000) return 5;
   if (x < 1000000) return 6;
   if (x < 10000000) return 7;
   if (x < 100000000) return 8;
   if (x < 1000000000) return 9;
   return 10;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Print a permutation
/// This function prints a permutation on the standard output using cycle notation. If @em name
/// is not 0, the name followed by an equal sign is printed before the permutation. For example,
/// the statement <tt>permPrint("Perm",P);</tt> could produce the following output:
/// <pre>
/// Perm=(1 9)(2 3 6)(4 5 7)
/// </pre>
/// Fixed points are always suppressed in the output.

void permPrint(const char *name, const Perm_t *perm)
{
   // check arguments
   permValidate(MTX_HERE, perm);

   // print the name
   if (name != NULL) {
      printf("%s=",name);
   }

   // print the permutation
   const uint32_t deg = perm->degree;
   const uint32_t* p = perm->data;

   uint8_t* done = NALLOC(uint8_t, deg);
   uint8_t* seed = done - 1;
   uint8_t* seedEnd = done + deg;

   // Run through all orbits. 
   int empty = 1;
   int count = 0;
   while (1) {
      // Find start of next orbit.
      while (*++seed != 0 && seed != seedEnd);
      if (seed == seedEnd)
         break;             // Done!
      uint32_t x = (uint32_t)(seed - done);

      // Suppress fixed points (GAP does not like them)
      if (p[x] == x) {
         done[x] = 1;
         continue;
      }

      empty = 0;
      int first = 1;
      while (1) {
         MTX_ASSERT_DEBUG(x >= 0 && x < deg);
         if (done[x])
            break;      // orbit complete
         done[x] = 1;

         if (first) {
            first = 0;
            if ((count += width(x) + 1) > 77) {
               printf("\n    (%lu",(unsigned long)x);
               count = 5 + width(x);
            } else {
               printf("(%lu",(unsigned long)x);
            }
         } else {
            if ((count += width(x)+1) > 77) {
               printf(",\n    %lu",(unsigned long)x);
               count = 4 + width(x);
            } else {
               printf(",%lu",(unsigned long)x);
            }
         }
         x = p[x];
      }
      printf(")");
      ++count;
   }

   sysFree(done);
   if (empty) {
      printf("()");
   }
   if (name != NULL) {
      printf("\n");
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Power of a permutation
/// This function calculates the n-th power of a permutation.
/// It allocates a new permutation, leaving the original
/// permutation intact. The caller is responsible for deleting the
/// result when it is no longer needed.
/// @param p Pointer to the permutation.
/// @param n Exponent. Must be greather than or equal to 0.
/// @return @em n-th power of @em p or 0 on error.

Perm_t *permPower(const Perm_t *p, int n)
{

   // check arguments
   permValidate(MTX_HERE, p);
   if (n < 0) {
      mtxAbort(MTX_HERE,"Invalid exponent %d < 0",n);
      return NULL;
   }

   Perm_t* q = permAlloc(p->degree);
   
   uint32_t* xp = p->data;
   uint32_t* xq = q->data;

   // calculate the n-th power
   for (uint32_t i = 0; i < p->degree; ++i) {
      uint32_t k = i;
      for (uint32_t l = n; l > 0; --l) {
         k = xp[k];
      }
      xq[i] = k;
   }
   return q;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void permConvertLegacyFormat(uint32_t *data, uint32_t degree)
{
   for (uint32_t i = 0; i < degree; ++i) {
      if (data[i] == 0) {
         return;
      }
   }
   for (int i = 0; i < degree; ++i) {
      --data[i];
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads permutation data from a file and returns the permutation.
/// Note: this function can only be called after a permutation header was read.
/// See also @ref permRead.
///
/// @param f File to read from.
/// @return Pointer to the permutation.

Perm_t* permReadData(MtxFile_t* f)
{
   const uint32_t objectType = mfObjectType(f);
   if (objectType != MTX_TYPE_PERMUTATION) {
      mtxAbort(MTX_HERE, "%s: bad type 0x%lx, expected 0x%lx (PERMUTATION)",
         f->name, (unsigned long) objectType, (unsigned long) MTX_TYPE_PERMUTATION);
   }

   Perm_t* p = permAlloc(f->header[1]);
   mfRead32(f, p->data, p->degree);
   permConvertLegacyFormat(p->data, p->degree);
   permValidate(MTX_HERE, p);

   // Make sure a second permReadData() call will fail.
   f->header[0] = 0xFFFFFFFF;

   return p;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads a permutation from a File.
/// This function reads a permutation from a file. The file must be opened for reading.
/// After return the file pointer is advanced to the first position after the permutation.
///
/// See also: @ref permLoad
///
/// @param f File to read from.
/// @return Pointer to the permutation, or 0 on error.

Perm_t* permRead(MtxFile_t* f)
{
   mfReadHeader(f);
   return permReadData(f);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads a single permutation from a file.
/// This function opens a file, reads a single permutation, closes the file, and returns the
/// permutation. If the file contains more than one permutation, only the first one is read.
/// 
/// See also: @ref permRead
///
/// @param fn File name.
/// @return Pointer to the permutation.

Perm_t* permLoad(const char* fn)
{
   MtxFile_t* f = mfOpen(fn, "rb");
   Perm_t* p = permRead(f);
   mfClose(f);
   return p;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes a permutation to a file. See also @ref permSave.

void permWrite(const Perm_t* perm, MtxFile_t* file)
{
   permValidate(MTX_HERE, perm);
   uint32_t hdr[3] = { MTX_TYPE_PERMUTATION, perm->degree, 1U };
   mfWrite32(file, hdr, 3);
   mfWrite32(file, perm->data, perm->degree);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes a permutation to a named file.
///
/// This function creates a file, writes a single permutation to the file and closes the file.
/// If a file with the same name already exists, its contents are destroyed.
/// See also @ref permWrite.

void permSave(const Perm_t* perm, const char* fileName)
{
   permValidate(MTX_HERE, perm);
   MtxFile_t* f = mfOpen(fileName, "wb");
   permWrite(perm, f);
   mfClose(f);
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
