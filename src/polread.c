////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Read a polynomial from a file
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


static long tmpfl = 0;
static long tmpdeg = 0;
static PTR tmpvec = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Allocate temporary workspace

static void mktmp(long fl, long deg)
{
   ffSetField(fl);
   if (deg > 0) {
   }
   if (tmpfl != fl || tmpdeg < deg) {
      if (tmpvec != NULL) { sysFree(tmpvec); }
      tmpvec = ffAlloc(1, deg + 1);
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
/// for deleting this polynomial with polFree() when it is no longer needed.
/// @see polLoad()
/// @param f File to read from.
/// @return Pointer to the polynomial, or 0 on error.

Poly_t *polRead(FILE *f)
{
   Poly_t *p;
   long hdr[3];

   if (sysReadLong32(f,hdr,3) != 3) {
      mtxAbort(MTX_HERE,"Cannot read header");
      return NULL;
   }
   if (hdr[0] != -2) {
      mtxAbort(MTX_HERE,"No polynomial (type=%d)",(int)hdr[0]);
      return NULL;
   }
   mktmp(hdr[1],hdr[2]);
   p = polAlloc(hdr[1],hdr[2]);
   if (p->Degree > 0) {
      int i;
      if (ffReadRows(f,tmpvec,1, p->Degree + 1) != 1) {
         polFree(p);
         mtxAbort(MTX_HERE,"Cannot read data");
         return NULL;
      }
      for (i = 0; i <= p->Degree; ++i) {
         p->Data[i] = ffExtract(tmpvec,i);
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
/// @see polRead()
/// @param fn File name.
/// @return Pointer to the polynomial read from the file, or 0 on error.

Poly_t *polLoad(const char *fn)
{
   FILE *f;
   Poly_t *p;

   if ((f = sysFopen(fn,"rb")) == NULL) {
      mtxAbort(MTX_HERE,"Cannot open %s",fn);
      return NULL;
   }
   p = polRead(f);
   fclose(f);
   if (p == NULL) {
      mtxAbort(MTX_HERE,"Cannot read polynomial from %s",fn);
      return NULL;
   }
   return p;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
