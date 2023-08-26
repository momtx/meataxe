////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Write a permutation
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


/// @addtogroup perm
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a Permutation to a File.
/// This function writes a permutation to a file. The file must be open
/// for writing.
/// @see PermSave
/// @param perm Permutation to write.
/// @param f File to write to.
/// @return The function returns 0 on success and -1 on error.

int permWrite(const Perm_t *perm, FILE *f)
{
   long hdr[3];

   permValidate(MTX_HERE, perm);
   hdr[0] = -1;
   hdr[1] = perm->Degree;
   hdr[2] = 1;
   if (sysWriteLong32(f,hdr,3) != 3) {
      mtxAbort(MTX_HERE,"Cannot write header");
      return -1;
   }
   if (sysWriteLong32(f,perm->Data,hdr[1]) != (int) hdr[1]) {
      mtxAbort(MTX_HERE,"Cannot write data");
      return -1;
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a Permutation to a File.
/// This function creates a file, writes a single permutation to the file and
/// closes the file. If a file with the specified name already exists, its
/// contents are destroyed.
/// @see PermWrite
/// @param perm Permutation to write.
/// @param fn File name.
/// @return The function returns 0 on success and -1 on error.

int permSave(const Perm_t *perm, const char *fn)
{
   FILE *f;
   int result;

   permValidate(MTX_HERE, perm);
   if ((f = sysFopen(fn,"wb")) == NULL) {
      mtxAbort(MTX_HERE,"Cannot open %s",fn);
      return -1;
   }
   result = permWrite(perm,f);
   fclose(f);
   if (result != 0) {
      mtxAbort(MTX_HERE,"Cannot write permutation to %s",fn);
      return -1;
   }
   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
