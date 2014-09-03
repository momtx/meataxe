/* ============================= C MeatAxe ==================================
   File:        $Id: ffio.c,v 1.2 2007-09-09 22:00:26 mringe Exp $
   Comment:     I/o for vectors and matrices.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include <string.h>
#include <stdlib.h>
#include "meataxe.h"

MTX_DEFINE_FILE_INFO

////////////////////////////////////////////////////////////////////////////////////////////////////

/// @addtogroup ff
/// @{

/// Read Rows
/// This function reads 1 or more rows from a file.
/// The row size must have been set before.
/// @param f Pointer to File.
/// @param buf Pointer to data buffer.
/// @param n Number of rows to read.
/// @return Number of rows that were actually read, or -1 on error.

int FfReadRows(FILE *f, PTR buf, int n)
{
   int i;
   register char *b = (char *) buf;

   /* Handle special case <FfNoc> = 0
      ------------------------------- */
   if (FfNoc == 0) {
      return n;
   }

   /* Read rows one by one
      -------------------- */
   for (i = 0; i < n; ++i) {
      if (fread(b,FfTrueRowSize(FfNoc),1,f) != 1) { break; }
      b += FfCurrentRowSize;
   }
   if (ferror(f)) {
      MTX_ERROR("Read failed: %S");
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
/// @return The number of rows that were successully written.
/// A return value different from |n| indicates an error.

int FfWriteRows(FILE *f, PTR buf, int n)
{
   int i;
   register char *b = (char *) buf;

   /* Handle special case <FfNoc> = 0
      ------------------------------- */
   if (FfNoc == 0) {
      return n;
   }

   /* Write rows one by one
      --------------------- */
   for (i = 0; i < n; ++i) {
      if (fwrite(b,FfTrueRowSize(FfNoc),1,f) != 1) { break; }
      b += FfCurrentRowSize;
   }
   if (ferror(f)) {
      MTX_ERROR("Write failed: %S");
   }
   return i;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Move to a Row.
/// This function sets the read/write pointer of file @em f to position
/// @em pos. I.e., the next FfReadRows() or FfWriteRows() will start
/// at the specified row.
/// Note that row numbers start with 0. If @em pos is different
/// from 0, the row size must have been set before with FfSetNoc().
///
/// You should always use FfSeekRow() rather than fseek() because
/// FfSeekRow() knows about MeatAxe file headers and adjusts the
/// file pointer appropriately.
/// @param f Pointer to File.
/// @param pos Row number (0-based).
/// @return 0 on success, -1 on error.

int FfSeekRow(FILE *f, int pos)
{
   long addr;

   if (FfOrder != -1) {
      addr = (long) FfTrueRowSize(FfNoc) * pos + 12;
   } else {
      addr = (long) FfNoc * 4 * pos + 12;
   }
   if (SysFseek(f,addr) == -1) {
      MTX_ERROR("Seek failed: %S");
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

FILE *FfReadHeader(const char *name, int *field, int *nor, int *noc)
{
   FILE *fd;
   long header[3];

   /* Open the file
      ------------- */
   fd = SysFopen(name,FM_READ);
   if (fd == NULL) {
      return NULL;
   }

   /* Read the file header
      -------------------- */
   if (SysReadLong(fd,header,3) != 3) {
      fclose(fd);
      MTX_ERROR1("%s: Error reading file header",name);
      return NULL;
   }

   /* Check header
      ------------ */
   if ((header[0] > 256) || (header[1] < 0) || (header[2] < 0)) {
      MTX_ERROR1("%s: Invalid header, possibly non-MeatAxe file",name);
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

FILE *FfWriteHeader(const char *name, int field, int nor, int noc)
{
   FILE *fd;
   long header[3];

   /* Open the file
      ------------- */
   fd = SysFopen(name,FM_CREATE);
   if (fd == NULL) {
      return NULL;
   }

   /* Write the file header
      --------------------- */
   header[0] = (long) field;
   header[1] = (long) nor;
   header[2] = (long) noc;
   if (SysWriteLong(fd,header,3) != 3) {
      MTX_ERROR1("%s: Error writing file header",name);
      fclose(fd);
      return NULL;
   }

   return fd;
}

/// @}
