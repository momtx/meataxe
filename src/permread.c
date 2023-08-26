////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Read a permutation
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup perm
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

void Perm_ConvertOld(long *data, int len)
{
   int i;

   // if point 0 exists the permutation is already in new format
   for (i = 0; i < len; ++i) {
      if (data[i] == 0) {
         return;
      }
   }

   /* Convert.
      --------  */
   for (i = 0; i < len; ++i) {
      --data[i];
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read a Permutation from a File.
/// This function reads a permutation from a file. @em f must be a pointer to an
/// open file with read permission.
/// If a permutation was successfully read, the function returns a pointer to
/// a newly created Perm_t object. The caller is responsible for deleting
/// this object as soon as the permutation is no longer needed.
/// @see permLoad()
/// @param f File to read from.
/// @return Pointer to the permutation, or 0 on error.

Perm_t *permRead(FILE *f)
{
   Perm_t *p;
   long hdr[3];

   if (sysReadLong32(f,hdr,3) != 3) {
      mtxAbort(MTX_HERE,"Cannot read header");
      return NULL;
   }
   if (hdr[0] != -1) {
      mtxAbort(MTX_HERE,"%s", MTX_ERR_NOTPERM);
      return NULL;
   }
   p = permAlloc(hdr[1]);
   if (p == NULL) {
      return NULL;
   }
   if (sysReadLong32(f,p->Data,p->Degree) != p->Degree) {
      permFree(p);
      mtxAbort(MTX_HERE,"Cannot read permutation data");
      return NULL;
   }
   Perm_ConvertOld(p->Data,p->Degree);
   permValidate(MTX_HERE, p);
   return p;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read a permutation.
/// This function opens a file, reads a single permutation, and closes the
/// file. The return value is a pointer to the permutation or 0 on
/// error. If the file contains more than one permutation, only the first one
/// is read.
/// If a permutation was successfully read, the function returns a pointer to
/// a newly created Perm_t object. The caller is responsible for deleting
/// this object as soon as the permutation is no longer needed.
/// @see permRead()
/// @param fn File name.
/// @return Pointer to the permutation read from the file, or 0 on error.

Perm_t *permLoad(const char *fn)
{
   FILE *f;
   Perm_t *p;

   if ((f = sysFopen(fn,"rb")) == NULL) {
      mtxAbort(MTX_HERE,"Cannot open %s",fn);
      return NULL;
   }
   p = permRead(f);
   fclose(f);
   if (p == NULL) {
      mtxAbort(MTX_HERE,"Cannot read permutation from %s",fn);
      return NULL;
   }
   return p;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
