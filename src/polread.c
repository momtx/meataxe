////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Read a polynomial from a file
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

MTX_DEFINE_FILE_INFO

static long tmpfl = 0;
static long tmpdeg = 0;
static PTR tmpvec = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Allocate temporary workspace

static void mktmp(long fl, long deg)
{
   FfSetField(fl);
   if (deg > 0) { FfSetNoc(deg + 1); }
   if ((tmpfl != fl) || (tmpdeg < deg)) {
      if (tmpvec != NULL) { SysFree(tmpvec); }
      tmpvec = FfAlloc(1);
      tmpdeg = deg;
      tmpfl = fl;
   }
}


/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read a polynomial from a file.
/// This function reads a polynomial from a file. If successful, the return
/// value is a pointer to a new Poly_t object. The caller is responsible
/// for deleting this polynomial with PolFree() when it is no longer needed.
/// @see PolLoad()
/// @param f File to read from.
/// @return Pointer to the polynomial, or 0 on error.

Poly_t *PolRead(FILE *f)
{
   Poly_t *p;
   long hdr[3];

   if (SysReadLong(f,hdr,3) != 3) {
      MTX_ERROR("Cannot read header");
      return NULL;
   }
   if (hdr[0] != -2) {
      MTX_ERROR1("No polynomial (type=%d)",(int)hdr[0]);
      return NULL;
   }
   mktmp(hdr[1],hdr[2]);
   if ((p = PolAlloc(hdr[1],hdr[2])) == NULL) {
      MTX_ERROR("Cannot create polynomial");
      return NULL;
   }
   if (p->Degree > 0) {
      int i;
      if (FfReadRows(f,tmpvec,1) != 1) {
         PolFree(p);
         MTX_ERROR("Cannot read data");
         return NULL;
      }
      for (i = 0; i <= p->Degree; ++i) {
         p->Data[i] = FfExtract(tmpvec,i);
      }
   }
   return p;
}


/// Read a Polynomial from a File.
/// This function opens a file, reads a single polynomial, and closes the
/// file. The return value is a pointer to the polynomial or NULL on
/// error. If the file contains more than one polynomial, only the first one
/// is read.
///
/// If a polynomial was successfully read, the function returns a pointer to
/// a newly created Poly_t object. The caller is responsible for deleting
/// this object as soon as it no longer needed.
/// @see PolRead()
/// @param fn File name.
/// @return Pointer to the polynomial read from the file, or 0 on error.

Poly_t *PolLoad(const char *fn)
{
   FILE *f;
   Poly_t *p;

   if ((f = SysFopen(fn,FM_READ)) == NULL) {
      MTX_ERROR1("Cannot open %s",fn);
      return NULL;
   }
   p = PolRead(f);
   fclose(f);
   if (p == NULL) {
      MTX_ERROR1("Cannot read polynomial from %s",fn);
      return NULL;
   }
   return p;
}


/// @}
