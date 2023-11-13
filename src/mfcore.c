////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Binary file I/O.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <errno.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


#define MF_MAGIC 0x229AE77B

/// @defgroup mf File I/O
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @class MtxFile_t
/// @brief
/// A MeatAxe binary file.
/// This structure serves as a handle for MeatAxe binary files with header and data part.

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Check a File Object.
/// This function checks if the argument points to a valid MtxFile_t structure.
/// @param file Pointer to the file.
/// @return 1 if @a file points to a valid file object, 0 otherwise.

int mfIsValid(const MtxFile_t *file)
{
   if (file == NULL)
      return 0;
   if (file->typeId != MF_MAGIC)
      return 0;
   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void mfValidate(const struct MtxSourceLocation* src, const MtxFile_t *file)
{
   if (file == NULL) {
      mtxAbort(src,"NULL file");
   }
   if (file->typeId != MF_MAGIC) {
      mtxAbort(src,"Invalid file");
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static MtxFile_t *allocFile(const char *name)
{
   MtxFile_t *f = ALLOC(MtxFile_t);
   memset(f,0,sizeof(*f));
   f->name = sysMalloc(strlen(name) + 1);
   strcpy(f->name,name);
   return f;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void freeFile(MtxFile_t *f)
{
   if (f->file != NULL) {
      fclose(f->file);
   }
   sysFree(f->name);
   memset(f,0,sizeof(*f));
   sysFree(f);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int isNonNegative(uint32_t x)
{
   return (x & 0x8000000L) == 0;
}

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
/// This function assumes that <tt>file->header</tt> contains a valid object header, i.e., that
/// @ref mfReadHeader or @ref mfTryReadHeader was called before.
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

/// Reads the object header at the current position.
/// On return, the members @c Field, @c Nor and @c Noc are set to the values found in the header.
/// The function fails and aborts the program if the header cannot be read or is invalid. This
/// includes the case that nothing could be read because the file pointer was already at the end of
/// file when the function was called.
/// See also @ref mfTryReadHeader.

void mfReadHeader(MtxFile_t* file)
{
   sysRead32(file->file, file->header, 3);
   checkHeader(file);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Tries to read an object header at the current position.
/// This function behaves like @ref mfReadHeader, but does not fail if the file pointer is already
/// at end of file.
/// @return 1 if an valid header was read, 0 otherwise (end of file or i/o error)

int mfTryReadHeader(MtxFile_t* file)
{
   if (!sysTryRead32(file->file, file->header, 3))
      return 0;
   checkHeader(file);
   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Opens a file for reading.
/// Fails if the file does not exist, or cannot be opened.

MtxFile_t *mfOpen(const char *name)
{
   MtxFile_t *f = allocFile(name);
   f->typeId = MF_MAGIC;
   f->file = sysFopen(name,"rb");
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

   /* Write the file header.
      ---------------------- */
   uint32_t header[3];
   header[0] = /*f->field = */ field;
   header[1] = /*f->nor = */nor;
   header[2] = /*f->noc = */noc;
   sysWrite32(f->file,header,3);

   f->typeId = MF_MAGIC;
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

/// Reads row vectors from a file.
/// @param f Pointer to the file.
/// @param buf Data buffer.
/// @param nor Number of rows to read.
/// @param noc Row size (number of columns).

void mfReadRows(MtxFile_t *f, PTR buf, uint32_t nor, uint32_t noc)
{
   mfValidate(MTX_HERE, f);

   // Handle empty data.
   if (nor == 0 || noc == 0)
      return;

   ffSetField(f->header[0]); // TODO: remove?
   const size_t rowSizeUsed = ffRowSizeUsed(noc);
   const size_t rowSize = ffRowSize(noc);

   // Read rows.
   uint8_t *b = (uint8_t *) buf;
   for (uint32_t i = nor; i > 0; --i) {
      if (fread(b, rowSizeUsed, 1, f->file) != 1) {
         mtxAbort(MTX_HERE,"%s: read error: %s", f->name, strerror(errno));
      }
      b += rowSize;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void mfSkip(MtxFile_t *file, size_t nBytes)
{
   if (fseek(file->file, (long)nBytes, SEEK_CUR) != 0)
      mtxAbort(MTX_HERE, "%s: seek error: %s", file->name, strerror(errno));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Write rows to a file.
/// See also @ref ffWriteRows.
///
/// @param f Pointer to the file.
/// @param buf Data buffer.
/// @param nor Number of rows to write.
/// @param noc Row size (number of columns).

void mfWriteRows(MtxFile_t *f, PTR buf, uint32_t nor, uint32_t noc)
{
   register char *b = (char *) buf;

   mfValidate(MTX_HERE, f);

   // Handle empty data
   if (nor == 0 || noc == 0) {
      return;
   }

   // Write rows one by one
   const size_t rowSize = ffRowSize(noc);
   const size_t rowSizeUsed = ffRowSizeUsed(noc);

   for (uint32_t i = nor; i > 0; --i) {
      if (fwrite(b,rowSizeUsed,1,f->file) != 1) {
         break;
      }
      b += rowSize;
   }

   if (ferror(f->file)) {
      mtxAbort(MTX_HERE,"%s: Write failed: %s",f->name, strerror(errno));
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes 32-bit integers to a file.
/// This function writes an array of 32-bit integers from @a buffer into a data file.
/// Each integer is written in LSB first format to the file. See @ref sysWrite32 for details.
/// @param f Pointer to the file.
/// @param buf Data buffer.
/// @param count Number of 32-bit integers to write.

void mfWrite32(MtxFile_t *f, const void *buf, int count)
{
   mfValidate(MTX_HERE, f);
   sysWrite32(f->file,buf,count);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads 32-bit integers from a file.
/// This function reads @a count 32-bit integers from a data file into a buffer.
/// Each integer is converted from file format (little-endian) into native format.
/// See @ref sysRead32 for details.

void mfRead32(MtxFile_t* f, void* buf, int count)
{
   mfValidate(MTX_HERE, f);
   sysRead32(f->file,buf,count);
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
