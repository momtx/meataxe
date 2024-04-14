////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Binary file I/O.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <errno.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data

/// @defgroup mf File I/O
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @class MtxFile_t
/// @brief
/// A MeatAxe binary file.
/// This structure serves as a handle for MeatAxe binary files with header and data part.

////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Verifies that the given file is valid, aborts the program otherwise.

void mfValidate(const struct MtxSourceLocation* src, const MtxFile_t *file)
{
   if (file == NULL) {
      mtxAbort(src,"NULL file");
   }
   if (file->typeId != MTX_TYPE_BINFILE) {
      mtxAbort(src,"Invalid file (bad typeId 0x%lx", (unsigned long) file->typeId);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static MtxFile_t *allocFile(const char *fileName)
{
   MTX_ASSERT(fileName != NULL);
   MtxFile_t *f = (MtxFile_t*) mmAlloc(MTX_TYPE_BINFILE, sizeof(MtxFile_t));
   f->name = sysMalloc(strlen(fileName) + 1);
   strcpy(f->name, fileName);
   return f;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void freeFile(MtxFile_t *f)
{
   if (f->file != NULL) {
      fclose(f->file);
      f->file = NULL;
   }
   sysFree(f->name);
   f->name = NULL;
   mmFree(f, MTX_TYPE_BINFILE);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int isNonNegative(uint32_t x)
{
   return (x & 0x8000000L) == 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int isValidFieldOrder(uint32_t x)
{
   #if MTX_ZZZ == 0
   static const uint32_t MTX_MAX_Q = 256;
   #elif MTX_ZZZ == 1
   static const uint32_t MTX_MAX_Q = 65536;
   #else
   #error
   #endif
   return x <= MTX_MAX_Q && x >= 2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int isValidHeader(uint32_t* type, const uint32_t header[3])
{
   if (!isNonNegative(header[1]) && isNonNegative(header[2])) {
      return 0;
   }
   *type = header[0];
   switch (header[0]) {
      case MTX_TYPE_PERMUTATION: return 1;
      case MTX_TYPE_POLYNOMIAL: return isValidFieldOrder(header[1]) && (int32_t) header[2] >= -1;
      case MTX_TYPE_BITSTRING_FIXED: return 1;
      case MTX_TYPE_BITSTRING_DYNAMIC: return 1;
      case MTX_TYPE_INTMATRIX: return 1;
      case MTX_TYPE_MATRIX: return 0; // never used, header[0] is field order
   }

   if (header[0] >= MTX_TYPE_BEGIN) {
      return 0;
   }

   *type = MTX_TYPE_MATRIX;
   return isValidFieldOrder(header[0]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the object type.
///
/// This function can only be called after @ref mfReadHeader returned or @ref mfTryReadHeader
/// returned zero. It returns the type id (e.g., MTX_TYPE_MATRIX) corresponding to the header that
/// has just been read.
///
/// The function fails and aborts the program if the header is invalid.

uint32_t mfObjectType(const MtxFile_t* file)
{
   uint32_t type = MTX_TYPE_BEGIN;
   if (!isValidHeader(&type, file->header))
      mtxAbort(MTX_HERE, "%s: invalid object header (0x%lx,0x%lx,0x%lx)",
            file->name,
            (unsigned long) file->header[0],
            (unsigned long) file->header[1],
            (unsigned long) file->header[2]);
   return type;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void checkHeader(MtxFile_t* f)
{
   mfObjectType(f);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads the object header at the current position and returns the object type.
///
/// The function fails and aborts the program if the header cannot be read or is invalid. This
/// includes the case that nothing could be read because the file pointer was already at the end of
/// file.
/// See also @ref mfTryReadHeader.

uint32_t mfReadHeader(MtxFile_t* file)
{
   sysRead32(file->file, file->header, 3);
   return mfObjectType(file);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Tries to read an object header at the current position.
///
/// This function behaves like @ref mfReadHeader, but does not fail if the file pointer is already
/// at end of file. The function can still fail if
/// - end of file was encountered after reading a part of the header, or
/// - any other i/o error.
///
/// @return 0 if an valid header was read, 1 otherwise (end of file)

int mfTryReadHeader(MtxFile_t* file)
{
   int result = mfTryRead32(file, file->header, 3);
   if (result == 1)
      return 1;
   checkHeader(file);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Opens a file.
/// Fails if the file does not exist or cannot be opened.
/// See @ref mfOpen for the meaning of @p mode.

MtxFile_t* mfOpen(const char* name, const char* mode)
{
   MtxFile_t* f = allocFile(name);
   f->file = sysFopen(name, mode);
   return f;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Opens a file for writing.
/// This functions creates a new file or truncates an existing file. The file is opened
/// for writing, and a MeatAxe file header is written to the file.

MtxFile_t *mfCreate(const char *name, uint32_t field, uint32_t nor, uint32_t noc)
{
   MtxFile_t *f = allocFile(name);
   f->file = sysFopen(name,"wb");

   uint32_t header[3];
   header[0] = field;
   header[1] = nor;
   header[2] = noc;
   sysWrite32(f->file,header,3);
   return f;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Closes a file and releases all associated memory.

void mfClose(MtxFile_t *file)
{
   mfValidate(MTX_HERE, file);
   freeFile(file);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Advances the file pointer by a given number of bytes.

void mfSkip(MtxFile_t *file, size_t nBytes)
{
   if (fseek(file->file, (long)nBytes, SEEK_CUR) != 0)
      mtxAbort(MTX_HERE, "%s: seek error: %s", file->name, strerror(errno));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes raw bytes to a file.
/// This function writes an array of 8-bit integers from a buffer to a data file.

void mfWrite8(MtxFile_t* f, const void* buf, size_t count)
{
   mfValidate(MTX_HERE, f);
   sysWrite8(f->file, buf, count);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads raw bytes from a file.
/// This function reads @p count bytes from a data file to a buffer.

void mfRead8(MtxFile_t* f, void* buf, size_t count)
{
   mfValidate(MTX_HERE, f);
   size_t nread = fread(buf, 1, count, f->file);
   if (nread != count) {
      const char *errorMessage = ferror(f->file) ? strerror(errno) : "unexpected end of file";
      mtxAbort(MTX_HERE, "Error reading %s: %s", f->name, errorMessage);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes 32-bit integers to a file.
/// This function writes an array of 32-bit integers from @p buf to a data file.
/// Each integer is written in LSB first format to the file. See @ref sysWrite32 for details.
/// @param f Pointer to the file.
/// @param buf Data buffer.
/// @param count Number of 32-bit integers to write.

void mfWrite32(MtxFile_t *f, const void *buf, size_t count)
{
   mfValidate(MTX_HERE, f);
   sysWrite32(f->file,buf,count);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads 32-bit integers from a file.
/// This function reads @p count 32-bit integers from a data file into a buffer.
/// Each integer is converted from file format (little-endian) into native format.
/// See @ref sysRead32 for details.

void mfRead32(MtxFile_t* f, void* buf, size_t count)
{
   mfValidate(MTX_HERE, f);
   // note: use mfTryRead32() instead of sysRead32() for better error messages
   int result = mfTryRead32(f, buf, count);
   if (result != 0) {
      mtxAbort(MTX_HERE, "Error reading %s: unexpected end of file", f->name);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads 32-bit integers from a file.
/// This function works like @ref mfRead32 but does not fail if the file pointer was already at the
/// end of file. The function still fails if ony a part of the requested data could be read before
/// the end of file was encountered, or if there was any other i/o error.
///
/// @return 0 on success (read complete), 1 on end-of-file (nothing as read)

int mfTryRead32(MtxFile_t* f, void* buf, size_t count)
{
   mfValidate(MTX_HERE, f);
   const int result = sysTryRead32(f->file, buf, count);
   if (result == 0 || result == 1) {
      return result;
   }
   mtxAbort(MTX_HERE, "Error reading %s: %s", f->name, strerror(errno));
   return -1;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
