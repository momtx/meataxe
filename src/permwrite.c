////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Write a permutation
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

/// @addtogroup perm
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes a permutation to a file. See also @ref permSave.

void permWrite(const Perm_t *perm, FILE *f)
{
   permValidate(MTX_HERE, perm);
   uint32_t hdr[3] = {MTX_TYPE_PERMUTATION, perm->Degree, 1U};
   sysWrite32(f,hdr,3);
   sysWrite32(f,perm->Data,perm->Degree);
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes a permutation to a file.
/// This function creates a file, writes a single permutation to the file and closes the file.
/// If a file with the same name already exists, its contents are destroyed.
/// See also @re3f PermWrite.

void permSave(const Perm_t *perm, const char *fileName)
{
   permValidate(MTX_HERE, perm);
   FILE* f = sysFopen(fileName,"wb");
   permWrite(perm,f);
   fclose(f);
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
