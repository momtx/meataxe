////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - I/o for vectors and matrices.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <string.h>
#include <stdlib.h>
#include "meataxe.h"



/// @addtogroup ff
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read Rows
/// This function reads 1 or more rows from a file.
/// The row size must have been set before.
/// @param f Pointer to File.
/// @param buf Pointer to data buffer.
/// @param n Number of rows to read.
/// @return Number of rows that were actually read, or -1 on error.

int ffReadRows(FILE *f, PTR buf, int n, int noc)
{
   int i;
   register char *b = (char *) buf;

   // handle special case: NOC=0
   if (noc == 0) {
      return n;
   }

   // read rows one by one
   for (i = 0; i < n; ++i) {
      if (fread(b,ffTrueRowSize(noc),1,f) != 1) {
	  break;
      }
      b += ffRowSize(noc);
   }
   if (ferror(f)) {
      mtxAbort(MTX_HERE,"Read failed: %S");
      return -1;
   }
   return i;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write rows
/// This function writes 1 or more rows to a file.
/// The row size must have been set before.
/// @param f Pointer to File.
/// @param buf Pointer to data buffer.
/// @param n Number of rows to write.
/// @return The number of rows (n), or -1 on error.

int ffWriteRows(FILE *f, PTR buf, int n, int noc)
{
   int i;
   register char *b = (char *) buf;

   // handle special case: NOC=0
   if (noc == 0) {
      return n;
   }

   // write rows one by one
   for (i = 0; i < n; ++i) {
      if (fwrite(b,ffTrueRowSize(noc),1,f) != 1) {
	  break;
      }
      b += ffRowSize(noc);
   }
   if (ferror(f)) {
      mtxAbort(MTX_HERE,"Write failed: %S");
      return -1;
   }
   return i;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Move to a Row.
/// This function sets the read/write pointer of file @em f to position
/// @em pos. I.e., the next ffReadRows() or ffWriteRows() will start
/// at the specified row.
/// Note that row numbers start with 0. If @em pos is different
/// from 0, the row size must have been set before with ffSetNoc().
///
/// You should always use ffSeekRow() rather than fseek() because
/// ffSeekRow() knows about MeatAxe file headers and adjusts the
/// file pointer appropriately.
/// @param f Pointer to File.
/// @param pos Row number (0-based).
/// @return 0 on success, -1 on error.

int ffSeekRow(FILE *f, int pos)
{
   long addr;

   if (ffOrder != -1) {
      addr = (long) ffTrueRowSize(ffNoc) * pos + 12;
   } else {
      addr = (long) ffNoc * 4 * pos + 12;
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
/// @param name File name.
/// @param field Pointer to buffer for first header entry (usually the field order).
/// @param nor Pointer to buffer for second header entry (usually the number of rows).
/// @param noc Pointer to buffer for third header entry (usually the number of columns).
/// @return Pointer to open file, or 0 on error.

FILE *ffReadHeader(const char *name, int *field, int *nor, int *noc)
{
   FILE *fd;
   long header[3];

   /* Open the file
      ------------- */
   fd = sysFopen(name, "rb");
   if (fd == NULL) {
      return NULL;
   }

   /* Read the file header
      -------------------- */
   if (sysReadLong32(fd,header,3) != 3) {
      fclose(fd);
      mtxAbort(MTX_HERE,"%s: Error reading file header",name);
      return NULL;
   }

   /* Check header
      ------------ */
#if MTXZZZ == 0
   #define MTX_MAX_Q 256
#elif MTXZZZ == 1
   #define MTX_MAX_Q 65535
#else
   #error
#endif
   if ((header[0] > MTX_MAX_Q) || (header[1] < 0) || (header[2] < 0)) {
      mtxAbort(MTX_HERE,"%s: Invalid header, possibly non-MeatAxe file",name);
      fclose(fd);
      return NULL;
   }

   /* Store header values in user-supplied buffers
      -------------------------------------------- */
   *field = (int)header[0];
   *nor = (int)header[1];
   *noc = (int)header[2];

   return fd;
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
/// @return Pointer to open file, or |NULL| on error.

FILE *ffWriteHeader(const char *name, int field, int nor, int noc)
{
   FILE *fd;
   long header[3];

   /* Open the file
      ------------- */
   fd = sysFopen(name,"wb");
   if (fd == NULL) {
      return NULL;
   }

   /* Write the file header
      --------------------- */
   header[0] = (long) field;
   header[1] = (long) nor;
   header[2] = (long) noc;
   if (sysWriteLong32(fd,header,3) != 3) {
      mtxAbort(MTX_HERE,"%s: Error writing file header",name);
      fclose(fd);
      return NULL;
   }

   return fd;
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
