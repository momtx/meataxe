////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Read a permutation
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup perm
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

void permConvertLegacyFormat(uint32_t *data, uint32_t degree)
{
   for (uint32_t i = 0; i < degree; ++i) {
      if (data[i] == 0) {
         return;
      }
   }
   for (int i = 0; i < degree; ++i) {
      --data[i];
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads a permutation using a given object header.
/// This function reads the permutation data (points), assuming the header has already been read
/// and is passed as second argument.
///
/// @param f File to read from.
/// @param header The object header.
/// @return Pointer to the permutation.

Perm_t *permReadData(FILE *f, const uint32_t header[3])
{
   if (header[0] != MTX_TYPE_PERMUTATION) {
      mtxAbort(MTX_HERE,"%s",MTX_ERR_NOTPERM);
   }
   Perm_t* p = permAlloc(header[1]);
   sysRead32(f, p->Data, p->Degree);
   permConvertLegacyFormat(p->Data, p->Degree);
   permValidate(MTX_HERE, p);
   return p;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Reads a permutation from a File.
/// This function reads a permutation from a file. The file must be opened for reading.
/// After return the file pointer is advanced to the firts position after the permtation.
///
/// See also: @ref permLoad
///
/// @param f File to read from.
/// @return Pointer to the permutation, or 0 on error.

Perm_t *permRead(FILE *f)
{
   uint32_t header[3];
   sysRead32(f,header,3);
   return permReadData(f, header);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads a single permutation from a file.
/// This function opens a file, reads a single permutation, closes the file, and returns the
/// permutation. If the file contains more than one permutation, only the first one is read.
/// 
/// See also: @ref permRead
///
/// @param fn File name.
/// @return Pointer to the permutation.

Perm_t *permLoad(const char *fn)
{
   int context = mtxBegin("Reading permutation: %s", fn);

   FILE *f;
   if ((f = sysFopen(fn,"rb")) == NULL) {
      mtxAbort(MTX_HERE,"Cannot open %s",fn);
   }
   Perm_t* p = permRead(f);
   fclose(f);

   mtxEnd(context);
   return p;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
