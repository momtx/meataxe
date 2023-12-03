////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - String formatting
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static const uint32_t TYPEID_STRING_BUFFER = 0x3B628F15LU;

/// @defgroup str String Builder
/// @{
/// These functions are used to work with dynamically allocated strings.

////////////////////////////////////////////////////////////////////////////////////////////////////

static int sbIsValid(const StrBuffer* sb)
{
   if (sb == NULL)
      return 0;
   if (sb->typeId != TYPEID_STRING_BUFFER)
      return 0;
   if (sb->size > sb->capacity)
      return 0;
   if (sb->capacity > 0 && sb->data == NULL)
      return 0;
   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void sbValidate(const struct MtxSourceLocation* sl, const StrBuffer* sb)
{
   if (!sbIsValid(sb))
      mtxAbort(sl, "Invalid string builder");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Allocates a new string builder.

StrBuffer* sbAlloc(size_t initialCapacity)
{
   StrBuffer* sb = ALLOC(StrBuffer);
   memset(sb, 0, sizeof(*sb));
   sb->capacity = initialCapacity;
   sb->data = NALLOC(char, initialCapacity + 1);
   *sb->data = 0;
   sb->typeId = TYPEID_STRING_BUFFER;
   return sb;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Deletes a string builder and releases all associated memory.

void sbFree(StrBuffer *sb)
{
   sbValidate(MTX_HERE, sb);
   sysFree(sb->data);
   memset(sb, 0, sizeof(*sb));
   sysFree(sb);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the string builder data.
/// The returned pointer becomes invalid if the string is modified, for example with @ref sbAppend.

const char* sbData(StrBuffer *sb)
{
   sbValidate(MTX_HERE, sb);
   char *copy = NALLOC(char, sb->size + 1);
   memcpy(copy, sb->data, sb->size + 1);
   return copy;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns a copy of the string builder data.
/// The returned pointer is dynamically allocated and remains valid after the string builder is
/// destroyed. The caller must release the string with sysFree() if it no longer needed.

char* sbCopy(StrBuffer* sb)
{
   sbValidate(MTX_HERE, sb);
   char* copy = NALLOC(char, sb->size + 1);
   memcpy(copy, sb->data, sb->size + 1);
   return copy;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Clears the string builder data. The string builder remains valid.

void sbClear(StrBuffer* sb)
{
   sbValidate(MTX_HERE, sb);
   sb->size = 0;
   MTX_ASSERT(sb->data != NULL);
   *sb->data = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the string builder data and destroys the string builder.
/// This is a (slightly more efficient) variant of sbCopy() + sbFree(). Like with @ref sbCopy(), the
/// caller becomes the owner the returned string and must release it with sysFree().

char* sbToString(StrBuffer* sb)
{
   char* data = sb->data;
   sb->data = NULL;
   sb->size = 0;
   sb->capacity = 0;
   sbFree(sb);
   return data;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns an ephemeral string containing the string builder data and destroys the string buffer.
/// The returned string is managed automatically and must not be release by the caller. 
/// See @ref strMakeEphemeral.

char* sbToEphemeralString(StrBuffer* sb)
{
   char* data = strMakeEphemeral(sb->data);
   sb->data = NULL;
   sb->size = 0;
   sb->capacity = 0;
   sbFree(sb);
   return data;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Makes sure that at least minFree characters (plus the terminating NUL) can be appended.

static void reserve(StrBuffer* sb, size_t minFree)
{
   if (sb->size + minFree > sb->capacity) {
      sb->capacity = sb->size + ((minFree < 100) ? 100 : minFree);
      sb->data = NREALLOC(sb->data, char, sb->capacity + 1);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Appends a fixed text fragment to the string.

void sbAppend(StrBuffer* sb, const char* fragment)
{
   sbValidate(MTX_HERE, sb);
   if (fragment == NULL)
      return;
   const size_t fragmentLength = strlen(fragment);
   reserve(sb, fragmentLength);
   memcpy(sb->data + sb->size, fragment, fragmentLength + 1);
   sb->size += fragmentLength;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Appends a formatted string with multiple arguments (vprintf style) to the string.

void sbVprintf(StrBuffer* sb, const char* fmt, va_list args)
{
   sbValidate(MTX_HERE, sb);
   reserve(sb, 100);
   va_list args2;
   va_copy(args2, args);
   int len = vsnprintf(sb->data + sb->size, sb->capacity - sb->size + 1, fmt, args2);
   va_end(args2);
   if (len < 0) {
      mtxAbort(MTX_HERE, "String formatting error");
   }
   if (len > 100) {
      reserve(sb, len);
      int len2 = vsnprintf(sb->data + sb->size, sb->capacity - sb->size + 1, fmt, args);
      MTX_ASSERT(len2 == len);
   }
   sb->size += len;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Appends a formatted string with multiple arguments (printf style) to the string.

void sbPrintf(StrBuffer* sb, const char* fmt, ...)
{
   va_list args;
   va_start(args, fmt);
   sbVprintf(sb, fmt, args);
   va_end(args);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Formats a string.
/// This function is similar to vsnprintf() but returns a pointer to a dynamically allocated buffer,
/// which must be relased by the caller.

char *strVMprintf(const char* s, va_list args)
{
   StrBuffer* sb = sbAlloc(100);
   sbVprintf(sb, s, args);
   return sbToString(sb);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Formats a string.
/// This function works like @ref strVMprintf but excepts a variable list of arguments instead of a
/// @c va_list.

char *strMprintf(const char* s, ...)
{
   va_list args;
   va_start(args, s);
   return strVMprintf(s, args);
   va_end(args);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Formats a string.
/// This function is similar to vsnprintf() but returns a pointer to an ephemeral string, which is
/// managed internally and will eventually released. See @ref strMakeEphemeral.

char *strVEprintf(const char* s, va_list args)
{
   StrBuffer* sb = sbAlloc(100);
   sbVprintf(sb, s, args);
   return sbToEphemeralString(sb);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Formats a string.
/// This function works like @ref strVEprintf but excepts a variable list of arguments instead of a
/// @c va_list.

char *strEprintf(const char* s, ...)
{
   va_list args;
   va_start(args, s);
   return strVEprintf(s, args);
   va_end(args);
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
