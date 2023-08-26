////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Error handling
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

const char MTX_ERR_NOMEM[] = "Not enough memory";
const char MTX_ERR_GAME_OVER[] = "Time limit exceeded";
const char MTX_ERR_DIV0[] = "Division by 0 or singular Matrix";
const char MTX_ERR_FILEFMT[] = "Bad format";
const char MTX_ERR_BADARG[] = "Bad argument";
const char MTX_ERR_RANGE[] = "Out of range";
const char MTX_ERR_NOTECH[] = "Matrix not in chelon form";
const char MTX_ERR_NOTSQUARE[] = "Matrix not square";
const char MTX_ERR_INCOMPAT[] = "Arguments are incompatible";
const char MTX_ERR_BADUSAGE[] = "Bad command line";
const char MTX_ERR_OPTION[] = "Bad usage of option";
const char MTX_ERR_NARGS[] = "Bad number of arguments";
const char MTX_ERR_NOTMATRIX[] = "Not a matrix";
const char MTX_ERR_NOTPERM[] = "Not a permutation";

////////////////////////////////////////////////////////////////////////////////////////////////////
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

typedef void ErrorHandler_t(const struct MtxErrorInfo *err);

static void (*ErrorHandler)(const struct MtxErrorInfo *err) = NULL;
static FILE *LogFile = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////

static void DefaultHandler(const struct MtxErrorInfo *e)
{
   if (LogFile == NULL) {
      LogFile = stderr;
   }

   if (e->source && e->source->file) {
      const char *baseName = strrchr(e->source->file, '/');
      baseName = baseName != NULL ? baseName + 1 : e->source->file;
      fprintf(LogFile,"%s:%d (%s): ", baseName, e->source->line, e->source->func);
   }
   fprintf(LogFile,"%s\n",e->message);
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
/// @param sl The source location to be included in the error message. Use the MTX_HERE to use
///    print the current location. If @a sl is NULL, no location will be shown.

void mtxAbort(const struct MtxSourceLocation* sl, const char *text, ...)
{
   va_list args;
   va_start(args,text);
   char txtbuf[2000];
   vsnprintf(txtbuf, sizeof(txtbuf), text, args);
   va_end(args);

   struct MtxErrorInfo err = { .source = sl, .message = txtbuf };
   if (ErrorHandler)
      ErrorHandler(&err);
   else
      DefaultHandler(&err);
   exit(9);
}

/// @}


// vim:fileencoding=utf8:sw=3:ts=8:et:cin
