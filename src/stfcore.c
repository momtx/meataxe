////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Structured Text File (STF) basic functions.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @defgroup stf Text File Handling
/// @{
/// The MeatAxe library provides functions for input and output of data in
/// human-readable text format. Files that are created with this set of functions
/// have a defined structure, and are refered to as "structured text files" (STF).
/// An examples for this type of files is the .cfinfo file which is
/// used by the submodule lattice programs.
///
/// <b>File format</b><br>
/// A structured text file is interpreted as a sequence of lines. While the STF input functions
/// can read very long lines, the output functions try to limit the line length to 80 characters
/// in order to make the file more readable.
/// Each line is one of the following:
/// - Lines starting with a "#" in column 1 are comment lines and are ignored completely.
///   Empty lines are ignored, too.
/// - A non-comment line with a non-blank character in column 1 marks the beginning of a new
///   entry. Such a line has the format
///   @code Name := Value@endcode
///   Both @em Name and @em Value are arbitrary strings, except that they cannot contain leading
///   or trailing blanks. In fact leading and trailing  blanks as well as any blanks around the
///   ":=" are removed on input.
/// - Lines starting with a whitspace character are interpreted as continuing lines.
///   Obviously a continuing line may occur only after an entry has started. The contents of
///   the continuing line, after leading lanks have been removed, are appended to @em Value.
///
/// <b>Data formats</b><br>
/// Besides the removal of leading and trailing blanks there is no restriction
/// on the format of <Value> in an STF entry. There are, however, predefined
/// functions that read and write integers and sequences of integers. An application
/// should use these functions where possible. The format used by the integer i/o
/// functions is most easily demonstrated in an example:
/// <pre>
/// Field := 7;
/// Multiplicity := [1,1,1,2];
/// Dimensions := [11,22,33,44,55];
/// </pre>
/// The format has been chosen such that GAP can read the text file without modification.

////////////////////////////////////////////////////////////////////////////////////////////////////

static void validateBase(const struct MtxSourceLocation* src, const StfData* f)
{
   if (f == NULL) {
      mtxAbort(src, "NULL file");
   }
   if (f->typeId != MTX_TYPE_STFILE) {
      mtxAbort(src, "Invalid text file (type=0x%lu)", (unsigned long) f->typeId);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Aborts the program if the passed STF object is not valid.

void stfValidate(const struct MtxSourceLocation* where, const StfData* f)
{
   validateBase(where, f);
   if (f->fileName == NULL) {
      mtxAbort(where, "Invalid text file (name=NULL)");
   }
   if (f->file == NULL) {
      mtxAbort(where, "Invalid text file (file=NULL)");
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Close a structured text file.
/// Closing the file implies that the memory occupied by the StfData structure is freed.
/// Thus, after return, @p f is invalid and must not be dereferenced.
/// @param f Pointer to an open structured text file (STF) object.

void stfClose(StfData *f)
{
   validateBase(MTX_HERE, f);
   mtxEnd(f->context);
   sysFree(f->fileName);
   if (f->file != NULL) {
      fclose(f->file);
      f->file = NULL;
   }
   if (f->lineBuf != NULL) {
      sysFree(f->lineBuf);
      f->lineBuf = NULL;
   }
   mmFree(f, MTX_TYPE_STFILE);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static const char* provideFilePosition(void* userData)
{
   StfData* f = (StfData*) userData;
   static char buffer[300];
   snprintf(buffer, sizeof(buffer), "at %s, line %d", f->fileName, f->lineNo);
   return buffer;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Opens a structured text file.
/// @p name and @p mode have the same semantics as with sysFopen().
///
/// Note: each call of @c stfOpen() creates a log context, which is closed in the corresponding
/// call of @ref stfClose. Applications creating log contexts must make sure that calls of
/// @ref mtxBegin / @ref mtxEnd call are properly nested with @ref stfOpen / @ref stfClose.

StfData *stfOpen(const char *name, const char* mode)
{
   StfData *f = (StfData*) mmAlloc(MTX_TYPE_STFILE, sizeof(StfData));
   f->lineBufSize = 250;
   f->lineBuf = NALLOC(char,f->lineBufSize);
   f->fileName = strdup(name);
   f->file = sysFopen(name, mode);
   f->context = mtxBeginP(provideFilePosition, f);
   return f;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
