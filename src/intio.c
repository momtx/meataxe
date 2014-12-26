/* ============================= C MeatAxe ==================================
   File:        $Id $
   Comment:     Input/output of integers.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include "config.h"

MTX_DEFINE_FILE_INFO

/// @addtogroup os
/// @{


/// Read long integers.
/// This function reads @ n long integers from the file @a f into the array 
/// @a buf. @a buf must point to a memory area of at least n*sizeof(long) 
/// bytes and the file must be open for reading. The return value indicates how
/// many integers have actually been read. This number may be less than
/// @a n because the end of file was encountered while reading. A negative
/// return value indicates a file i/o error.
/// 
/// %SysReadLong32() expects that the numbers in the file are 4-byte integers 
/// in little-endian format, i.e. the least significant byte first.
/// Using a machine-independent data format makes MeatAxe data files 
/// more portable, but there are also some disadvantages:
/// - The conversion to and from machine-independent format involves several
///   arithmetic operations for each number read/written.
/// - The highest number which can be read/written is 2<sup>32</sup>-1.
/// @param f File to read from.
/// @param buf Pointer to buffer.
/// @param n Number of integers to read.
/// @return Number of integers that were actually read, or $/1$ on error.

int SysReadLong32(FILE *f, long *buf, int n)
{
    unsigned char a[4];
    int nread;

    /* Check the arguments
       ------------------- */
    if (f == NULL || buf == NULL || n < 0)
    {
	MTX_ERROR3("Invalid arguments (f:%s,buf:%s,n=%d)",f ? "ok" : "NULL",
		buf ? "ok" : "NULL",n);
	return -1;
    }

    /* Read <n> long integers
       ---------------------- */
    if (Mtx_IsX86)
    {
	nread = (int)fread((char *)buf,sizeof(long),(size_t)n,f);
	/* printf("nread=%d, error=%d, eof=%d, GetLastError()=%d\n",
	   nread,ferror(f),feof(f),GetLastError());
	*/
    }
    else
    {
    	for (nread = 0; nread < n; ++nread)
    	{
	    if (fread(a,1,4,f) != 4) break;
	    *buf = ((unsigned long)a[0]|((unsigned long)a[1] << 8)|
		((unsigned long)a[2] << 16)|((long)(char)a[3] << 24));
	    ++buf;
    	}
    }

    /* Check for errors
       ---------------- */
    if (ferror(f) && !feof(f))
    {
	MTX_ERROR("Read failed: %S");
	return -1;
    }
    return nread;
}



/// Write long integers.
/// This function writes @a n long integers from the the array @a buf to the 
/// file @a f. @a buf must point to a memory area of at least n*sizeof(long) 
/// bytes and @a f must be open for writing. The numbers are written in a 
/// machine-independent format which can be read by SysReadLong().
/// @param f File to write to.
/// @param buf Pointer to buffer.
/// @param n Number of integers to write.
/// @return Number of integers that were written, or $-1$ on error.

int SysWriteLong32(FILE *f, const long *buf, int n)
{
    unsigned char a[4];
    int nwritten;

    if (Mtx_IsX86)
	nwritten = fwrite((char *)buf,sizeof(long),n,f);
    else
    {
    	for (nwritten = 0; nwritten < n; ++nwritten)
    	{
    	    a[0] = (unsigned char) *buf;
    	    a[1] = (unsigned char) (*buf >> 8);
    	    a[2] = (unsigned char) (*buf >> 16);
    	    a[3] = (unsigned char) (*buf >> 24);
    	    if (fwrite(a,1,4,f) != 4) break;
	    ++buf;
    	}
    }
    return nwritten;
}


#if MTX_CONFIG_BIG_ENDIAN
static void Swap(long *dest, const long *src, int n)
{
    for (; n > 0; --n)
    {
	register long x = *src++;
#if MTX_CONFIG_LONG32
	*dest++ =
	      (x << 24)
	    + ((x << 8) & 0x00FF0000L)
	    + ((x >> 8) & 0x0000FF00L)
	    + (x >> 24);
#elif MTX_CONFIG_LONG64
	*dest++ =
	      (x << 56)
	    + ((x << 40) & 0x00FF000000000000L)
	    + ((x << 24) & 0x0000FF0000000000L)
	    + ((x <<  8) & 0x000000FF00000000L)
	    + ((x >>  8) & 0x00000000FF000000L)
	    + ((x >> 24) & 0x0000000000FF0000L)
	    + ((x >> 40) & 0x000000000000FF00L)
	    + (x >> 56);
#else
#error "sizeof(long) not known - check config.h"
#endif
    }
}
#endif

#define BLK_SIZE 4096

int SysReadLongX(FILE *f, long *buf, int n_bytes)
{
    int n_read = 0;		/* Bytes read so far */
#if MTX_CONFIG_BIG_ENDIAN
    unsigned long tmp[BLK_SIZE];
    while (n_read < num_elem)
    {
	int n_blk = sizeof(tmp);
	if (n_read + n_blk > n_bytes)
	    n_blk = n_bytes - n_read;
	if (fread(tmp,1,n_blk,f) != n_blk)
	    break;
	Swap(buf + n_read / sizeof(long),tmp,(n_blk - 1)/sizeof(long) + 1);
	n_written += n_blk;
    }
#else
    n_read = fread(buf,1,n_bytes,f);
#endif
    return n_read;
}


int SysWriteLongX(FILE *f, const long *buf, int n_bytes)
{
    int n_written = 0;
#if MTX_CONFIG_BIG_ENDIAN
    unsigned long tmp[BLK_SIZE];
    while (n_written < n_bytes)
    {
	int n_blk = sizeof(tmp);
	if (n_written + n_blk > n_bytes)
	    n_blk = n_bytes - n_written;
	Swap(tmp,buf + n_written / sizeof(long),(n_blk - 1) / sizeof(long) + 1);
	if (fwrite(tmp,1,n_blk,f) != n_blk)
	    break;
	n_written += n_blk;
    }
#else
    n_written = fwrite(buf,1,n_bytes,f);
#endif
    return n_written;
}



/// @}

