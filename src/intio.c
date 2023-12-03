////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Input/output of integers
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <string.h>
#include <errno.h>


/// @addtogroup os
/// @{

/// Reads an array of 32-bit integers.
/// This function expects that integers are stored in little-endian format in the file.
///
/// @param f Input file.
/// @param buf Pointer to data buffer (uint32_t* or int32_t*).
/// @param n Number of 32-bit integers to read.

void sysRead32(FILE *f, void* buf, size_t n)
{
   MTX_ASSERT(n >= 0);
   size_t nRead = fread(buf, 4, n, f);
   if (nRead != n)
      mtxAbort(MTX_HERE, "File read error: %s", strerror(errno));
   if (mtxIsBigEndian()) {
      uint8_t* p = (uint8_t*) buf;
      uint8_t* end = p + 4 * n;
      while (p < end) {
         uint8_t tmp = p[0];
         p[0] = p[3];
         p[3] = tmp;
         tmp = p[1];
         p[1] = p[2];
         p[2] = tmp;
         p += 4;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Tries to read an array of 32-bit integers
/// Returns 0 on success (read complete), 1 on end-of-file (nothing read), -1 on partial read
/// and -2 on any other error.

int sysTryRead32(FILE *f, void* buf, size_t n)
{
   MTX_ASSERT(n >= 0);
   size_t nRead = fread(buf, 4, n, f);
   if (nRead == 0 && feof(f)) {
      // TODO: this may be wrong if there were 1..3 bytes available 
      // TODO: clarify the exact behaviour of fread if the first item is incomplete
      // TODO: mayby use (buf,1,4*n,f) instead to avoid this problem
      return 1;
   }
   if (nRead != n)
      return ferror(f) ? -2 : -1;
     
   if (mtxIsBigEndian()) {
      uint8_t* p = (uint8_t*) buf;
      uint8_t* end = p + 4 * n;
      while (p < end) {
         uint8_t tmp = p[0];
         p[0] = p[3];
         p[3] = tmp;
         tmp = p[1];
         p[1] = p[2];
         p[2] = tmp;
         p += 4;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes an array of 32-bit integers.
/// This function expects that integers are stored in little-endian format in the file.
///
/// @param f Output file.
/// @param buf Pointer to data buffer (uint32_t* or int32_t*).
/// @param n Number of 32-bit integers to read.

void sysWrite32(FILE *f, const void* buf, size_t n)
{
   MTX_ASSERT(n >= 0);
   uint8_t* temp = NULL;
   const uint8_t* bp = buf;

   if (mtxIsBigEndian()) {
      uint32_t* temp = NALLOC(uint32_t, n);
      const uint8_t* p = (const uint8_t*) buf;
      const uint8_t* end = p + 4 * n;
      uint8_t* t = (uint8_t*) temp;
      while (p < end) {
         *t++ = p[3];
         *t++ = p[2];
         *t++ = p[1];
         *t++ = p[0];
         p += 4;
         bp = (const uint8_t*) temp;
      }
   }
   while (n > 0) {
      size_t nWritten = fwrite(bp, 4, n, f);
      if (nWritten == 0)
         mtxAbort(MTX_HERE, "File write error: %s", strerror(errno));
      n -= nWritten;
   }
   if (mtxIsBigEndian())
      sysFree(temp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads an array of 16-bit integers.
/// This function expects that integers are stored in little-endian format in the file.
///
/// @param f Input file.
/// @param buf Pointer to data buffer (uint16_t* or int16_t*).
/// @param n Number of 16-bit integers to read.

void sysRead16(FILE *f, void* buf, size_t n)
{
   MTX_ASSERT(f != NULL);
   size_t nRead = fread(buf, 2, n, f);
   if (nRead != n)
      mtxAbort(MTX_HERE, "File read error: %s", strerror(errno));
   if (mtxIsBigEndian()) {
      uint8_t* p = (uint8_t*) buf;
      uint8_t* end = p + 2 * n;
      while (p < end) {
         uint8_t tmp = p[0];
         p[0] = p[1];
         p[1] = tmp;
         p += 2;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes an array of 16-bit integers.
/// This function expects that integers are stored in little-endian format in the file.
///
/// @param f Output file.
/// @param buf Pointer to data buffer (uint16_t* or int16_t*).
/// @param n Number of 16-bit integers to read.

void sysWrite16(FILE *f, const void* buf, size_t n)
{
   MTX_ASSERT(f != NULL);
   uint8_t* temp = NULL;
   const uint8_t* bp = buf;

   if (mtxIsBigEndian()) {
      uint16_t* temp = NALLOC(uint16_t, n);
      const uint8_t* p = (const uint8_t*) buf;
      const uint8_t* end = p + (size_t) 2 * n;
      uint8_t* t = (uint8_t*) temp;
      while (p < end) {
         *t++ = p[1];
         *t++ = p[0];
         p += 2;
         bp = (const uint8_t*) temp;
      }
   }
   while (n > 0) {
      size_t nWritten = fwrite(bp, 2, n, f);
      if (nWritten == 0)
         mtxAbort(MTX_HERE, "File write error: %s", strerror(errno));
      n -= nWritten;
   }
   if (mtxIsBigEndian())
      sysFree(temp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads bytes from a file. Aborts the program on error.
///
/// @param f Input file.
/// @param buf Pointer to data buffer.
/// @param n Number of bytes to read.

void sysRead8(FILE *f, void* buf, size_t n)
{
   MTX_ASSERT(n >= 0);
   size_t nRead = fread(buf, 1, n, f);
   if (nRead != n)
      mtxAbort(MTX_HERE, "File read error: %s", strerror(errno));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes bytes to a file. Aborts the program on error.
///
/// @param f Output file.
/// @param buf Pointer to data buffer.
/// @param n Number of bytes to read.

void sysWrite8(FILE *f, const void* buf, size_t n)
{
   MTX_ASSERT(n >= 0);
   uint8_t* temp = NULL;
   const uint8_t* bp = buf;

   while (n > 0) {
      size_t nWritten = fwrite(bp, 1, n, f);
      if (nWritten == 0)
         mtxAbort(MTX_HERE, "File write error: %s", strerror(errno));
      n -= nWritten;
   }
   if (mtxIsBigEndian())
      sysFree(temp);
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
