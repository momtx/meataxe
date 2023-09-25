////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Write polynomial into a file
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

/// @addtogroup poly
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes a polynomial to a file. See also @ref polSave.

void polWrite(const Poly_t *p, FILE *file)
{
   polValidate(MTX_HERE, p);

   uint32_t hdr[3] = { MTX_TYPE_POLYNOMIAL, p->field, p->degree };
   sysWrite32(file,hdr,3);
   ffSetField(p->field);
   if (p->degree >= 0) {
      PTR tmp = ffAlloc(1, p->degree + 1);
      for (unsigned i = 0; i <= p->degree; ++i)
         ffInsert(tmp, i, p->data[i]);
      ffWriteRows(file,tmp, 1, p->degree + 1);
      ffFree(tmp);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes a single polynomial to a file.
/// If a file with e same name alread eyists, its contents are destroyed.
/// See also @ref PolWrite
///
/// @param pol Polynomial to write.
/// @param fn File name.

void polSave(const Poly_t *pol, const char *fn)
{
   polValidate(MTX_HERE, pol);

   FILE *file = sysFopen(fn,"wb");
   polWrite(pol,file);
   fclose(file);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads a polynomial using a given object header.
///
/// This function reads the polynomial coefficients, assuming the header has already been read
/// and is passed as second argument.
///
/// @param f File to read from.
/// @param header The object header.
/// @return Pointer to the polynomial-

Poly_t *polReadData(FILE *f, const uint32_t header[3])
{
   if (header[0] != MTX_TYPE_POLYNOMIAL)
      mtxAbort(MTX_HERE,"Bad object type 0x%lx, expected polynomial",(unsigned long)header[0]);

   const int field = (int) header[1];
   const size_t degree = header[2];
   ffSetField(field);
   Poly_t* polynomial = polAlloc(field, degree);
   if (polynomial->degree >= 0) {
      PTR tmpvec = ffAlloc(1, degree + 1);
      ffReadRows(f, tmpvec, 1, degree + 1);
      for (size_t i = 0; i <= degree; ++i) {
         polynomial->data[i] = ffExtract(tmpvec,i);
      }
      ffFree(tmpvec);
   }
   return polynomial;
}

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
   uint32_t hdr[3];
   sysRead32(f,hdr,3);
   if (hdr[0] != MTX_TYPE_POLYNOMIAL) {
      mtxAbort(MTX_HERE,"No polynomial (type=%d)",(int)hdr[0]);
   }

   const int field = hdr[1];
   const int degree = hdr[2];
   ffSetField(field);
   Poly_t* p = polAlloc(field, degree);
   if (p->degree > 0) {
      PTR tmpvec = ffAlloc(1, degree + 1);
      ffReadRows(f, tmpvec, 1, degree + 1);
      for (size_t i = 0; i <= degree; ++i) {
         p->data[i] = ffExtract(tmpvec,i);
      }
      ffFree(tmpvec);
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

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
