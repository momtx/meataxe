////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Matrix representations, file i/o
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup mrep
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Save a Matrix Representation.
/// This function saves all generators of a matrix representation.
/// Each generator ist written to different file. The file name
/// is constructed by appending ".1", ".2" etc. to @a basename or, if
/// @a basename contains a "%d" placeholder, by replacing the "%d"
/// with "1", "2", etc.
/// @param rep Pointer to the matrix representation.
/// @param basename Base file name for generators.
/// @return 0 on success, -1 on error.

int mrSave(const MatRep_t *rep, const char *basename)
{
   char *fn;
   int ext_format;          /* '%d' found in <basename> */
   int i;

   /* Make a copy of the basename an reserve extra bytes for the extension
      -------------------------------------------------------------------- */
   fn = sysMalloc(strlen(basename) + 10);
   if (fn == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate buffer");
      return -1;
   }

   /* Write the generators.
      --------------------- */
   ext_format = strstr(basename,"%d") != NULL;
   for (i = 0; i < rep->NGen; ++i) {
      if (ext_format) {
         sprintf(fn,basename,i + 1);
      } else {
         sprintf(fn,"%s.%d",basename,i + 1);
      }
      matSave(rep->Gen[i],fn);
   }

   /* Clean up.
      --------- */
   sysFree(fn);
   return i >= rep->NGen ? 0 : -1;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
