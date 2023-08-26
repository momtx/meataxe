////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Write polynomial into a file
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
   if (deg > 0) { ffSetNoc(deg + 1); }
   if ((tmpfl != fl) || (tmpdeg < deg)) {
      if (tmpvec != NULL) { sysFree(tmpvec); }
      tmpvec = ffAlloc(1, deg + 1);
      tmpdeg = deg;
      tmpfl = fl;
   }
}


/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a polynomial to a file.
/// @see polSave()
/// @param p Pointer to the polynomial.
/// @param f File to write to.
/// @return 0 on success, -1 on error.

int polWrite(const Poly_t *p, FILE *f)
{
   long hdr[3];

   polValidate(MTX_HERE, p);
   mktmp(p->Field,p->Degree);
   hdr[0] = -2;
   hdr[1] = p->Field;
   hdr[2] = p->Degree;
   if (sysWriteLong32(f,hdr,3) != 3) {
      mtxAbort(MTX_HERE,"Cannot write header");
   }
   if (p->Degree >= 0) {
      int i;
      for (i = 0; i <= p->Degree; ++i) {
         ffInsert(tmpvec,i,p->Data[i]);
      }
      if (ffWriteRows(f,tmpvec,1, p->Degree + 1) != 1) {
         mtxAbort(MTX_HERE,"Cannot write data");
      }
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a polynomial to a File.
/// This function creates a file, writes a single polynomial to the file and
/// closes the file. If a f ile with the specified name already exists, it's
/// contents are destroyed.
/// @see PolWrite
/// @param pol Polynomial to write.
/// @param fn File name.
/// @return 0 on success, -1 on error.

int polSave(const Poly_t *pol, const char *fn)
{
   FILE *f;
   int result;

   polValidate(MTX_HERE, pol);
   if ((f = sysFopen(fn,"wb")) == NULL) {
      mtxAbort(MTX_HERE,"Cannot open %s",fn);
      return -1;
   }
   result = polWrite(pol,f);
   fclose(f);
   if (result != 0) {
      mtxAbort(MTX_HERE,"Cannot write polynomial to %s",fn);
      return -1;
   }
   return result;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
