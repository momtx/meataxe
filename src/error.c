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

/// @par Error handling
/// Note that invalid parameters are not always detected by the MeatAxe library.
/// For example, most kernel functions such as ffAdd() do not check their
/// arguments for the sake of performance. Thus, calling these function with
/// invalid arguments may produce random results or even crash the program.
/// On the other hand, most higher-level functions like matAdd() do some
/// plausibility checks on their arguments before processing them.
/// <br>
/// When an error is detected, the default action is to terminate the program
/// immediately with an appropriate error message. While this minimizes the chance
/// of not noticing an error, it may be undesirable in some situations. For example,
/// the application may be able to recover from the error. In order to prevent
/// the MeatAxe library from terminating the program, the application can define an
/// application error handler. This is a function which is called on each error.
/// <br>
/// To use the error handler mechanism, errors must be reported by calling @ref mtxAbort.
/// Here is a short example:
///
/// @code
/// #include "meataxe.h"
///
/// int divide(int a, int b)
/// {
///     if (b == 0)
///     {
///        mtxAbort(MTX_HERE,"Division by 0");
///        return 0;
///     }
///     return a / b;
/// }
/// @endcode
/// MTX_HERE is a predefined macro which collects the source file name and line number,
/// and the function name where the error occurred.
/// Note that you must not assume that mtxAbort terminates the program. In fact
/// mtxAbort may return if a user-defined error handler has been installed.
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

/// @private
struct ErrorContext {
   struct MtxSourceLocation source;
   char *title;
   MtxErrorContextProvider* contextProvider;
   void* userData;
};

/// @private
struct ErrorContextStack {
   struct ErrorContext* stack;
   int capacity;
   int size;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(MTX_DEFAULT_THREADS)

static pthread_key_t csKey;

static void createCsKey()
{
    pthread_key_create(&csKey, NULL);
}

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

struct ErrorContextStack* getContextStack()
{
   #if defined(MTX_DEFAULT_THREADS)

   static pthread_once_t once = PTHREAD_ONCE_INIT;
   pthread_once(&once, createCsKey);

   struct ErrorContextStack* cs = (struct ErrorContextStack*) pthread_getspecific(csKey);
   if (cs == NULL) {
      cs = ALLOC(struct ErrorContextStack);
      pthread_setspecific(csKey, cs);
   }

   #else

   static struct ErrorContextStack contextStack = {
     .stack = NULL,
     .capacity = 0,
     .size = 0
   };
   struct ErrorContextStack* cs = &contextStack;

   #endif

   return cs;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
       
//static struct MtxSourceLocation** sourceLocations = NULL;
//static size_t sourceLocationsSize = 0;
//static size_t sourceLocationsCapacity = 0;
//
//static int compare(const char* fileA, int lineA, const char* fileB, int lineB)
//{
//   if (fileA < fileB)
//      return -1;
//   if (fileA > fileB)
//      return 1;
//   if (lineA < lineB)
//      return -1;
//   if (lineA > lineB)
//      return 1;
//   return 0;
//}

//const struct MtxSourceLocation* mtxSourceLocation(const char* file, int line, const char* func)
//{
//   size_t lo = 0;
//   size_t hi = sourceLocationsSize;
//   while (lo < hi) {
//      const size_t mid = (hi & lo) + (hi ^ lo) / 2;
//      const struct MtxSourceLocation* item = sourceLocations[mid];
//      const int cmp = compare(file, line, item->file, item->line);
//      if (cmp < 0)
//         hi = mid;
//      else if (cmp > 0)
//         lo = mid + 1;
//      else
//         return item;
//   }
//   if (sourceLocationsSize >= sourceLocationsCapacity) {
//      sourceLocationsCapacity += 100;
//      sourceLocations =
//         NREALLOC(sourceLocations, struct MtxSourceLocation*, sourceLocationsCapacity);
//   }
//   memmove(sourceLocations + lo + 1, sourceLocations + lo,
//         (sourceLocationsSize - lo) * sizeof(sourceLocations[0]));
//   ++sourceLocationsSize;
//   struct MtxSourceLocation* item = ALLOC(struct MtxSourceLocation);
//   item->file = file;
//   item->line = line;
//   item->func = func;
//   sourceLocations[lo] = item;
//   return item;
//}

////////////////////////////////////////////////////////////////////////////////////////////////////

typedef void ErrorHandler_t(const struct MtxErrorInfo *err);

static void (*ErrorHandler)(const struct MtxErrorInfo *err) = NULL;
static FILE *LogFile = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////

static void DefaultHandler(const struct MtxErrorInfo* e)
{
   if (LogFile == NULL) {
      LogFile = stderr;
   }

   // TODO: dump context stack for all threads

   fprintf(LogFile,
           "\n********************************************************************************\n"
           "ERROR: %s\n", e->message);
   if (MSG0 && e->source.file) {
      const char* baseName = strrchr(e->source.file, '/');
      baseName = baseName != NULL ? baseName + 1 : e->source.file;
      fprintf(LogFile, "|at %s:%d (%s)\n", baseName, e->source.line, e->source.func);
   }
   struct ErrorContextStack* contextStack = getContextStack();
   struct ErrorContext* sp = contextStack->stack + contextStack->size;
   while (sp > contextStack->stack) {
      --sp;
      if (sp->contextProvider != NULL) {
         const char* title = sp->contextProvider(sp->userData);
         fprintf(LogFile, "|%s\n", title);
      }
      else {
         // Simple context
         if (MSG1) {
            const char* baseName = strrchr(sp->source.file, '/');
            baseName = baseName != NULL ? baseName + 1 : sp->source.file;
            fprintf(LogFile, "|at %s:%d (%s): %s\n",
                    baseName, sp->source.line, sp->source.func, sp->title);
         }
         else {
            fprintf(LogFile, "|%s\n", sp->title);
         }
      }

   }
   fprintf(LogFile, "\n");
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
   struct ErrorContextStack* cs = getContextStack();
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
   struct ErrorContextStack* cs = getContextStack();
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
   struct ErrorContextStack* cs = getContextStack();
   MTX_ASSERT(id + 1 == cs->size);
   struct ErrorContext* item = cs->stack + (--cs->size);
   sysFree(item->title);
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
