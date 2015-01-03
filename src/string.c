////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - String functions.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>
#include <stdarg.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
//   Local data

MTX_DEFINE_FILE_INFO
static size_t const EMPTY[4] = {0,0,0,0};

/// @defgroup str Dynamic Strings
/// @{
/// These functions are used to work with dynamically allocated strings.
/// A dynamic string contains a normal «char*» pointing to a NUL terminated
/// text.
/// Note however, that dynamic strings use their own memory management
/// which cannot be mixed with the standard C library memory functions.
/// Unused strings must be freed with StrFree(), and you must never use
/// free() or realloc() on a dynamic string.

#define REP(s) ((size_t *)(s) - 3)
#define LEN(s) (((size_t *)(s))[-1])
#define CAP(s) (((size_t *)(s))[-2])
#define REF(s) (((size_t *)(s))[-3])
#define TXT(h) ((char*)(h) + 3 * sizeof(size_t))
#define AVL(s) (CAP(s) - LEN(s))

////////////////////////////////////////////////////////////////////////////////////////////////////

static void safe_strcpy(char *d, const char *s, size_t len)
{
   (void) Mtx_ThisFile;
   if (len > 0) {
      if (s) { memcpy(d,s,len); }
      d[len] = 0;
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static char *mkrep(const char *s, size_t len, size_t cap)
{
   MTX_ASSERT(len <= cap);
   if ((len == 0) && (cap == 0)) { return TXT(EMPTY); }
   size_t *rep = (size_t *) SysMalloc(cap + 1 + 3 * sizeof(size_t));
   rep[0] = 1;
   rep[1] = cap;
   rep[2] = len;
   char *txt = TXT(rep);
   safe_strcpy(txt,s,len);
   return txt;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void delref(char *s)
{
   size_t * const r = REP(s);
   if ((r[0] != 0) && (--r[0] == 0)) {
      SysFree(r);
   }
}


/*
   static void addref(char *s)
   {
    size_t * const r = REP(s);
    if (r[0] > 0)
   ++r[0];
   }
 */

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Makes sure that: REF == 1, CAP >= «cap», CAP >= LEN

static void excl(char **s, unsigned cap)
{
   if (REF(*s) != 1) {                          /* Kopie erzeugen? */
      if (cap < LEN(*s)) { cap = LEN(*s); }
      char * const old_s = *s;
      *s = mkrep(*s,LEN(*s),cap);
      delref(old_s);
   } else if (cap > CAP(*s)) {                  /* Vergrößern? */
      size_t *rep = REP(*s);
      rep = (size_t *) SysRealloc(rep,3 * sizeof(size_t) + cap + 1);
      rep[1] = cap;
      *s = TXT(rep);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Makes at least «min_avail» bytes available without chaging LEN

static void reserve(char **s, size_t min_avail)
{
   excl(s,LEN(*s) + min_avail);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Create a string.
/// This function creates an empty string.
/// @see StrFree()
/// @param initial_capacity Number of bytes to reserve on allocation.
/// @return The new string.

String StrAlloc(size_t initial_capacity)
{
   String s;
   s.S = (initial_capacity > 0) ? mkrep(0,0,initial_capacity) : TXT(EMPTY);
   return s;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Free a string.
/// This function frees a dynamic string.
/// @see StrAlloc()

void StrFree(String *s)
{
   delref(s->S);
   s->S = 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void append(String *s, const char *src, size_t len)
{
   MTX_ASSERT(s->S != 0);
   size_t const my_len = LEN(s->S);
   if ((src >= s->S) && (src <= s->S + my_len)) {
      /* Note: excl() may reallocate the buffer, adjust src accordingly. */
      size_t const n = src - s->S;
      excl(&s->S,my_len + len);
      src = s->S + n;
   } else {
      excl(&s->S,my_len + len);
   }

   if (s->S != TXT(EMPTY)) {
      safe_strcpy(s->S + my_len,src,len);
      LEN(s->S) = my_len + len;
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Appends text to the end of a string.
/// @param s The string to be modified.
/// @param text Text to append.

void StrAppend(String *s, const char *text)
{
   return append(s, text, strlen(text));
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// vsnprintf() replacement.
/// Formats a string with vsnprintf() and appends the resulting
/// text to «s».
/// @param s The string to be modified.
/// @param fmt Format
/// @param args Arguments for the % placeholders in fmt.

void StrVAppendF(String *s, const char *fmt, va_list args)
{
   va_list saved_args;
   va_copy(saved_args,args);
   reserve(&s->S,strlen(fmt) + 32);
   int n = vsnprintf(s->S + LEN(s->S),AVL(s->S) + 1,fmt,args);
   if (n < 0) {
      return;
   }
   if ((size_t) n > AVL(s->S)) {
      reserve(&s->S,n);
      vsnprintf(s->S + LEN(s->S),n + 1,fmt,saved_args);
   }
   size_t const new_len = LEN(s->S) + n;
   LEN(s->S) = new_len;
   s->S[new_len] = 0;
   va_end(saved_args);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// snprintf() replacement.
/// Works like StrVAppendF() but expects a variable argument list.

void StrAppendF(String *s, const char *fmt, ...)
{
   va_list al;
   va_start(al,fmt);
   StrVAppendF(s,fmt,al);
   va_end(al);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// vsnprintf() replacement.
/// Works like StrVAppendF() but overwrites the string with the formatted text.
/// @param s The string to be modified.
/// @param fmt Format
/// @param args Arguments for the % placeholders in fmt.

void StrVPrintF(String *s, const char *fmt, va_list args)
{
   va_list saved_args;
   va_copy(saved_args,args);
   reserve(&s->S, strlen(fmt) + 32);
   int n = vsnprintf(s->S,CAP(s->S) + 1,fmt,args);
   if (n < 0) {
      return;
   }
   if ((size_t) n > CAP(s->S)) {
      reserve(&s->S, n - LEN(s->S));
      vsnprintf(s->S,n + 1,fmt,saved_args);
   }
   LEN(s->S) = n;
   va_end(saved_args);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// snprintf() replacement.
/// Works like StrVPrintF() but expects a variable argument list.

void StrPrintF(String *s, const char *fmt, ...)
{
   va_list al;
   va_start(al,fmt);
   StrVPrintF(s,fmt,al);
   va_end(al);
}


/// @}