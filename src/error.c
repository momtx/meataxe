////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Error handling
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/// @addtogroup app
/// @{

const char MTX_ERR_GAME_OVER[] = "Time limit exceeded";
const char MTX_ERR_DIV0[] = "Division by 0 or singular Matrix";
const char MTX_ERR_FILEFMT[] = "Bad format";
const char MTX_ERR_BADARG[] = "Bad argument";
const char MTX_ERR_RANGE[] = "Out of range";
const char MTX_ERR_NOTECH[] = "Matrix not in chelon form";
const char MTX_ERR_NOTSQUARE[] = "Matrix not square";
const char MTX_ERR_INCOMPAT[] = "Arguments are incompatible";
const char MTX_ERR_OPTION[] = "Bad usage of option";
const char MTX_ERR_NOTMATRIX[] = "Not a matrix";
const char MTX_ERR_NOTPERM[] = "Not a permutation";

////////////////////////////////////////////////////////////////////////////////////////////////////
       
typedef void ErrorHandler_t(const struct MtxErrorInfo *err);
static void (*ErrorHandler)(const struct MtxErrorInfo *err) = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////

static void DefaultHandler(const struct MtxErrorInfo* e)
{
   logPrepareForAbort();

   // TODO: dump context stack for all threads
   logPrintf(MTX_LOG_ERROR, "**********************************************************");
   logPrintf(MTX_LOG_ERROR, "FATAL ERROR: %s", e->message);
   if (e->source.file) {
      const char* baseName = strrchr(e->source.file, '/');
      baseName = baseName != NULL ? baseName + 1 : e->source.file;
      logPrintf(MTX_LOG_ERROR, "|at %s:%d (%s)", baseName, e->source.line, e->source.func);
   }
   struct ErrorContextStack* contextStack = pexContextStack();
   struct ErrorContext* sp = contextStack->stack + contextStack->size;
   while (sp > contextStack->stack) {
      --sp;
      if (sp->contextProvider != NULL) {
         const char* title = sp->contextProvider(sp->userData);
         logPrintf(MTX_LOG_ERROR, "|%s", title);
      }
      else {
         const char* baseName = strrchr(sp->source.file, '/');
         baseName = baseName != NULL ? baseName + 1 : sp->source.file;
         logPrintf(MTX_LOG_ERROR, "|at %s:%d (%s): %s",
               baseName, sp->source.line, sp->source.func, sp->title);
      }
   }
   exit(9);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Define an application error handler.
/// This function defines an application error handler which is called every
/// time an error occurs inside the MeatAxe library.
/// @param h Address of the new error handler or NULL to restore the default handler.
/// @return The previous error handler.

MtxErrorHandler_t *MtxSetErrorHandler(MtxErrorHandler_t *h)
{
   MtxErrorHandler_t *old = ErrorHandler;
   ErrorHandler = h;
   return old;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Terminates the program with an error message. 
/// @param sl The source location to be included in the error message. Pass @c MTX_HERE to use
///    print the current location. If @a sl is NULL, no location will be shown.
/// @param text The error message (printf style format), followed by any arguments.

void mtxAbort(const struct MtxSourceLocation* sl, const char *text, ...)
{
   va_list args;
   va_start(args,text);
   char txtbuf[2000];
   vsnprintf(txtbuf, sizeof(txtbuf), text, args);
   va_end(args);

   struct MtxErrorInfo err = { .message = txtbuf };
   if (sl != NULL) {
      err.source = *sl;
   } else {
      memset(&err.source, 0, sizeof(err.source));
   }
   if (ErrorHandler)
      ErrorHandler(&err);
   DefaultHandler(&err);
   exit(9);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Adds error context information.

int mtxBegin(const struct MtxSourceLocation* sl, const char *s, ...)
{
   struct ErrorContextStack* cs = pexContextStack();
   if (cs->size >= cs->capacity) {
      cs->capacity = cs->size + 10;
      cs->stack = NREALLOC(cs->stack, struct ErrorContext, cs->capacity);
   }
   struct ErrorContext* item = cs->stack + cs->size;
   memset(item, 0, sizeof(*item));
   va_list args;
   va_start(args, s);
   item->source = *sl;
   item->title = strVMprintf(s, args);
   va_end(args);
   return cs->size++;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int mtxBeginScope(MtxErrorContextProvider ec, void* userData)
{
   struct ErrorContextStack* cs = pexContextStack();
   if (cs->size >= cs->capacity) {
      cs->capacity = cs->size + 10;
      cs->stack = NREALLOC(cs->stack, struct ErrorContext, cs->capacity);
   }
   struct ErrorContext* item = cs->stack + cs->size;
   memset(item, 0, sizeof(*item));
   //item->source = ...
   item->title = NULL;
   item->contextProvider = ec;
   item->userData = userData;
   return cs->size++;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Removes error context information.

void mtxEnd(int id)
{
   struct ErrorContextStack* cs = pexContextStack();
   MTX_ASSERT(id + 1 == cs->size);
   struct ErrorContext* item = cs->stack + (--cs->size);
   sysFree(item->title);
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
