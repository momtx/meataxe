////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Write a permutation
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

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

int PermWrite(const Perm_t *perm, FILE *f)
{
   long hdr[3];

   if (!PermIsValid(perm)) {
      return -1;
   }
   hdr[0] = -1;
   hdr[1] = perm->Degree;
   hdr[2] = 1;
   if (SysWriteLong(f,hdr,3) != 3) {
      MTX_ERROR("Cannot write header");
      return -1;
   }
   if (SysWriteLong(f,perm->Data,hdr[1]) != (int) hdr[1]) {
      MTX_ERROR("Cannot write data");
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

int PermSave(const Perm_t *perm, const char *fn)
{
   FILE *f;
   int result;

   if (!PermIsValid(perm)) {
      return -1;
   }
   if ((f = SysFopen(fn,FM_CREATE)) == NULL) {
      MTX_ERROR1("Cannot open %s",fn);
      return -1;
   }
   result = PermWrite(perm,f);
   fclose(f);
   if (result != 0) {
      MTX_ERROR1("Cannot write permutation to %s",fn);
      return -1;
   }
   return 0;
}


/// @}
