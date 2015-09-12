////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Read a permutation
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
/// @see PermLoad()
/// @param f File to read from.
/// @return Pointer to the permutation, or 0 on error.

Perm_t *PermRead(FILE *f)
{
   Perm_t *p;
   long hdr[3];

   if (SysReadLong(f,hdr,3) != 3) {
      MTX_ERROR("Cannot read header");
      return NULL;
   }
   if (hdr[0] != -1) {
      MTX_ERROR1("%E", MTX_ERR_NOTPERM);
      return NULL;
   }
   p = PermAlloc(hdr[1]);
   if (p == NULL) {
      return NULL;
   }
   if (SysReadLong(f,p->Data,p->Degree) != p->Degree) {
      PermFree(p);
      MTX_ERROR("Cannot read permutation data");
      return NULL;
   }
   Perm_ConvertOld(p->Data,p->Degree);
   if (!PermIsValid(p)) {
      return NULL;
   }
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
/// @see PermRead()
/// @param fn File name.
/// @return Pointer to the permutation read from the file, or 0 on error.

Perm_t *PermLoad(const char *fn)
{
   FILE *f;
   Perm_t *p;

   if ((f = SysFopen(fn,FM_READ)) == NULL) {
      MTX_ERROR1("Cannot open %s",fn);
      return NULL;
   }
   p = PermRead(f);
   fclose(f);
   if (p == NULL) {
      MTX_ERROR1("Cannot read permutation from %s",fn);
      return NULL;
   }
   return p;
}


/// @}
