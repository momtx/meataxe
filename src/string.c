////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - String utilities
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Check if a string starts with a given prefix.
/// If @p s starts with @p prefix, the function returns a pointer to the caracter in @p s that 
/// follows the prefix. If @p s and @p prefix are equal, the return value points to the terminating
/// NUL byte of @p s.
/// If @p s does not start with @p prefix, the function returns NULL.
/// If either pointer is NULL, the function returns NULL

const char* strPrefix(const char* s, const char* prefix)
{
   if (s == NULL || prefix == NULL)
       return NULL;
   const size_t prefixLen = strlen(prefix);
   return (strncmp(s, prefix, prefixLen) == 0) ? s + prefixLen : NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns a copy of a character sequence.
/// The returned pointer is a NUL terminated string containing a copy of the range
/// [@p begin, @p end). If either argument is NULL or if @p begin > @p end, the function fails.

char* strRange(const char *begin, const char* end)
{
   MTX_ASSERT(begin != NULL);
   MTX_ASSERT(end != NULL);
   MTX_ASSERT(end >= begin);
   size_t n = end - begin;
   char *range = NALLOC(char, n + 1);
   memcpy(range, begin, end - begin);
   range[n] = 0;
   return range;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns an independent copy of a string.

char* strDup(const char *src)
{
    MTX_ASSERT(src != NULL);
    return strRange(src, src + strlen(src));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Compares two character sequences lexicographically.
/// Both sequences are given by their begin and end pointers. The end pointer may be NULL meaning
/// that the sequence extends from begin up to the first NUL character.
/// If any begin pointer is NULL or begin > end, the function fails.

int strCompareRange(const char* x, const char* xEnd, const char* y, const char* yEnd)
{
   MTX_ASSERT(x != NULL);
   MTX_ASSERT(y != NULL);
   if (xEnd == NULL) { xEnd = x + strlen(x); } else { MTX_ASSERT(x <= xEnd); }
   if (yEnd == NULL) { yEnd = y + strlen(y); } else { MTX_ASSERT(y <= yEnd); }

   while (1) {
      if (x == xEnd) {
         return (y == yEnd) ? 0 : -1;
      }
      if (y == yEnd) {
         return 1;
      }
      if (*x < *y) { return -1; }
      if (*x > *y) { return 1; }
      ++x;
      ++y;
   }
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
