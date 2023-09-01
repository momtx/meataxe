////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - I/o for vectors and matrices.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <errno.h>
#include <stdlib.h>
#include <string.h>

/// @addtogroup ff
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read Rows
/// This function reads 1 or more rows from a file.
/// @param f Pointer to File.
/// @param buf Pointer to data buffer.
/// @param nor Number of rows to read.
/// @param noc Row size (number of columns).

void ffReadRows(FILE *f, PTR buf, int nor, int noc)
{
   if (noc == 0 || nor == 0) {
      return;
   }
   const size_t rowSizeUsed = ffRowSizeUsed(noc);
   const size_t rowSize = ffRowSize(noc);

   // read rows one by one
   char *b = (char *) buf;
   for (int i = 0; i < nor; ++i) {
      if (fread(b,rowSizeUsed,1,f) != 1) {
         mtxAbort(MTX_HERE,"Error reading row(s) from file: %s", strerror(errno));
      }
      b += rowSize;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write rows
/// This function writes 1 or more rows to a file.
/// The row size must have been set before.
/// @param f Pointer to File.
/// @param buf Pointer to data buffer.
/// @param nor Number of rows to write.
/// @param noc Row size (number of columns).

void ffWriteRows(FILE *f, PTR buf, int nor, int noc)
{
   if (noc == 0 || nor == 0) {
      return;
   }
   size_t rowSizeUsed = ffRowSizeUsed(noc);
   size_t rowSize = ffRowSize(noc);


   const char *b = (const char *) buf;
   for (int i = 0; i < nor; ++i) {
      if (fwrite(b,rowSizeUsed,1,f) != 1) {
         mtxAbort(MTX_HERE,"Error writing row(s) to file: %s", strerror(errno));
      }
      b += rowSize;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Move to a Row.
/// This function sets the read/write pointer of file @em f to position
/// @em pos. I.e., the next ffReadRows() or ffWriteRows() will start
/// at the specified row.
/// Note that row numbers start with 0.
///
/// You should always use ffSeekRow() rather than fseek() because
/// ffSeekRow() knows about MeatAxe file headers and adjusts the
/// file pointer appropriately.
/// @param f Pointer to File.
/// @param pos Row number (0-based).
/// @param noc Row size (number of columns).
/// @return 0 on success, -1 on error.

int ffSeekRow(FILE *f, int pos, int noc)
{
   long addr;

   if (ffOrder != -1) {
      addr = ffRowSizeUsed(noc) * pos + 12;
   } else {
      addr = (long) noc * 4 * pos + 12;
   }
   if (sysFseek(f,addr) == -1) {
      mtxAbort(MTX_HERE,"Seek failed: %S");
      return -1;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Open File and Read Header.
/// This function opens a data file for input and reads the file header (3 integers).
/// The exact meaning of the header values depends on the file type.
/// For a matrix they are field order, number of rows and number of columns.
/// See @ref pg_file_formats for details.
/// @param fileName File name.
/// @param field Pointer to buffer for first header entry (usually the field order).
/// @param nor Pointer to buffer for second header entry (usually the number of rows).
/// @param noc Pointer to buffer for third header entry (usually the number of columns).
/// @return Pointer to open file.

FILE *ffReadHeader(const char *fileName, int *field, int *nor, int *noc)
{
   FILE *file = sysFopen(fileName, "rb");
   int32_t header[3];
   sysRead32(file, header, 3);

   /* Check header
      ------------ */
#if MTX_ZZZ == 0
   #define MTX_MAX_Q 256
#elif MTX_ZZZ == 1
   #define MTX_MAX_Q 65535
#else
   #error
#endif
   // TODO: check header[1..2] for all known types
   if (header[0] > MTX_MAX_Q) {
      mtxAbort(MTX_HERE,"%s: Invalid header, possibly not a MeatAxe file",fileName);
      fclose(file);
      return NULL;
   }

   /* Store header values in user-supplied buffers
      -------------------------------------------- */
   *field = (int)header[0];
   *nor = (int)header[1];
   *noc = (int)header[2];

   return file;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Open file and write header.
/// This function opens a data file for input and writes the the file header.
/// If the file does not exist, a new file is created. If the file exists it
/// is overwritten.
/// @param name File name.
/// @param field First header entry (usually the field order).
/// @param nor Second header entry (usually the number of rows).
/// @param noc Third header entry (usually the number of columns).
/// @return Pointer to open file.

FILE *ffWriteHeader(const char *name, int field, int nor, int noc)
{
   FILE *fd = sysFopen(name,"wb");
   uint32_t header[3] = {(uint32_t) field, (uint32_t) nor, (uint32_t) noc};
   sysWrite32(fd,header,3);
   return fd;
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
