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
/// Load a Matrix Representation.
/// This function creates a new matrix representation and reads the generators
/// from files. Each generator ist expected in a different file. The file name
/// is constructed by appending ".1", ".2" etc. to @a basename or, if @a basename
/// contains a "%d" placeholder, by replacing the "%d" with "1", "2", etc.
/// For example, the following lines
/// @code
/// m11 = mrLoad("m11",2);
/// m11 = mrLoad("m11.%d",2);
/// @endcode
/// are equivalent. In both cases, two matrices are read from "m11.2" and "m11.2",
/// repectively.
/// @param basename Base file name for generators.
/// @param ngen Number of generators.
/// @return Pointer to the representation.

MatRep_t *mrLoad(const char *basename, int ngen)
{
   char *fn;
   int ext_format;          /* '%d' found in <basename> */
   MatRep_t *mr;
   int i;

   /* Make a copy of the basename and reserve extra bytes for the extension.
      ---------------------------------------------------------------------- */
   fn = sysMalloc(strlen(basename) + 10);
   if (fn == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate buffer");
   }

   mr = mrAlloc(0,NULL,0);

   /* Read the generators
      ------------------- */
   ext_format = strstr(basename,"%d") != NULL;
   for (i = 0; i < ngen; ++i) {
      Matrix_t *gen;
      if (ext_format) {
         sprintf(fn,basename,i + 1);
      } else {
         sprintf(fn,"%s.%d",basename,i + 1);
      }
      gen = matLoad(fn);
      mrAddGenerator(mr,gen,0);
   }

   sysFree(fn);
   return mr;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
