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
/// Close a structured text file.
/// This function closes a structured text file.
/// Closing the file implies that the memory occupied by the StfData structure is freed.
/// Thus, after return, @a f is invalid and must not be dereferenced.
/// @param f Pointer to an open structured text file (STF) object.
/// @return 0 on success, non-zero on error.

int stfClose(StfData *f)
{
   if (f == NULL) {
      return -1;
   }
   if (f->File != NULL) {
      fclose(f->File);
   }
   if (f->LineBuf != NULL) {
      free(f->LineBuf);
   }
   memset(f,0,sizeof(StfData));
   free(f);
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Initialize a STF data structure
/// This funcion initializes a StfData structure.
///  The function assures that either the data structure is completely
///  initialized or, in case of any error, no resources are allocated.
/// @return 0 on success, -1 on error

static int stfInitData(StfData *f)
{
   memset(f,0,sizeof(StfData));
   f->LineBufSize = 250;
   f->LineBuf = NALLOC(char,f->LineBufSize);
   if (f->LineBuf == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate line buffer");
      return -1;
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Open a structured text file.
/// This function opens a structured text file. It returns to a StfData structure which
/// can be used with stfXxx() functions.
/// @a name and @a mode have the same semantics as with sysFopen().
///
/// @return StfData structure or 0 on error.

StfData *stfOpen(const char *name, const char* mode)
{
   StfData *f;

   // allocate and initialize a new STF data structure
   f = ALLOC(StfData);
   if (f == NULL) {
      return NULL;
   }
   if (stfInitData(f) != 0) {
      free(f);
      return NULL;
   }

   // open the file
   f->File = sysFopen(name, mode);
   if (f->File == NULL) {
      stfClose(f);
      return NULL;
   }

   return f;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
