////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Finite field arithmetic and common functions
// This is the "Small" version for field orders (q <= 256). Originally based on hprout.c written
// by Klaus Lux.
//
// (C) Copyright 1998-2014 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
// Modifications by Tom Hoffman, Mathematics Department, University of Arizona.
// Contributions by Simon King <simon.king@uni-jena.de>
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#include "meataxe.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Gobal data


typedef unsigned char BYTE;
static int MPB = 0;             /* No. of marks per byte */

int mtx_subfields[17];         // public list of proper subfields, terminated with 0

static int isFel(FEL x) { return (unsigned int) x < (unsigned int) ffOrder; }


/// @defgroup ff Finite Fields
///
/// The finite field part of the kernel provides finite field arithmetic and
/// basic operations with vectors and matrices over finite fields.
/// The kernel cannot operate simultaneously with different finite fields,
/// because there is a global row size and a global field order which must
/// be maintained by the higher layers.
///
/// There are two finite field modules available: one for small fields (up
/// to 256) and one for larger fields (up to 2<sup>16</sup>). The finite field
/// module is selected at compile time.
///
/// @par 'Small' Kernel (q≤256)
/// In the "small" kernel, field elements of GF(q) are represented by the numbers 0,1,...,q-1.
/// The field is defined by its Conway polynomial p(x)∈ℤₚ[x], where q=pⁿ (see [@ref HEBA05]).
/// Thus, we have a one-to-one correspondence of field elements a∈GF(q)
/// and polynomials fₐ(x)∈ℤₚ[x] of degree less than n.
/// By treating fₐ(x) as a polynomial over ℤ we can calculate f_a(p), giving the number of
/// the field element a. It follows that the elements of the prime field are represented by
/// 0,1,...p-1. The number 0 represents the zero element, and 1 represents the unit element.
/// <br>
/// For small fields (q <= 16), rows of matrices are compressed by packing two or more elements
/// in a byte. For example, if q=5, three elements (a,b,c) are combined into a single byte as
/// a + 5*b + 25*c.
///
/// @par 'Big' Kernel (q≤65536)
/// The big version stores field elements in 16-bit integers, i.e., each field element
/// occupies two bytes. Non-zero elements are stored as their logarithms with
/// respect to a fixed generator. In particular, the unit element is represented by
/// the integer 0. The zero element is represented by the special value 0xFFFF.
///
/// As a consequence of the different representations of field elements
/// in the small and big version, there are some rules which should be
/// respected by all programs:
/// - Never assign numbers to variables of type FEL or pass
///   numbers to functions expecting an argument of type FEL.
/// - Never perform integer arithmetic on variables of type FEL.
/// - Never printf() or scanf() variables of type FEL.
/// - Never use the literals 0 or 1 where the zero and unit element
///   of the field is intended. Instead you should use the FF_ZERO
///   and FF_ONE constants define in |meataxe.h|.
/// - Do not cast an integer to FEL or vice versa. Use
///   ffFromInt() and ffToInt() instead.
///
/// @par FEL and PTR
/// The kernel defines two basic data types:
/// FEL represents a single field element, and PTR is a pointer to a row vector.
/// PTR may be defined as <tt>*FEL</tt>, but this is  not mandatory.
/// <br>
/// The kernel also defines two constants:
/// FF_ZERO is the zero element of the current field, and
/// FF_ONE is unit element of the current field.
/// Depending on which kernel you are using, FF_ZERO and FF_ONE may be constants,
/// or they may be defined as variables or function calls.
///
/// @par Embedding of subfields
/// The field F=GF(p^n) has, for each m dividing n, a (unique) subfield G<F with p^m
/// elements. However, the numeric representation of field elements in F and G is different.
/// For example, the elements of the subfield of order 9 in GF(81) are
/// 0, 1, 2, 73, 74, 72, 38, 36, and 37, while the elements of GF(9) are represented as
/// 0, 1, ..., 8.
/// <br>
/// There are two functions, @ref ffEmbed and @ref ffRestrict, which convert FEL values
/// between a field and its subfields. The following example code converts a vector over
/// GF(9) to GF(27). Since the MeatAxe cannot handle two fields at the same time, it is
/// necessary to unpack the row over GF(9), change to GF(81) and pack the embedded elements
/// into a new row vector.
/// @code
/// ffSetField(9);
/// PTR row9 = ffAlloc(10);
/// // [...] - fill the row with values
/// FEL buf[10];
/// for (int i = 0; i &lt; 10; ++i) buf[i] = ffExtract(row9,i);
/// ffSetField(81);
/// PTR row81 = ffAlloc(10);
/// for (i = 0; i &lt; 10; ++i) ffInsert(row81,i,ffEmbed(buf[i],9));
/// @endcode

/// @addtogroup ff
/// @par Row vectors and row operations
/// A (row) vector is stored in memory as an array of bytes. The actual
/// size of a row depends on the number of marks in the row and on the
/// field order, because several field elements may be packed into one
/// byte. However, the size is always a multiple of sizeof(long).
/// Unused bytes at the end of a row (and even the unused bits in partially
/// used bytes) must always be filled with zeroes.
/// Otherwise, some functions as  ffCmpRows() may produce wrong results.
/// <br>
/// Packing of field elements is used for field orders less than 17. Let q be
/// the field order and m the largest natural number with q<sup>m</sup>≤256.
/// Then, m elements k<sub>0</sub>,...k<sub>n-1</sub> are packed into one byte
/// as k<sub>0</sub>+k<sub>1</sub>q+k<sub>2</sub>q<sup>2</sup>+...
/// Packing of field elements is used exclusively for rows (vectors),
/// not for polynomials or other data types.
/// <br>
/// Because of packing you cannot access the marks of a row with the usual bracket notation.
/// Instead, the MeatAxe library provides functions which store and extract marks.
/// <br>
/// Packed rows must be initialized before they are used. ffInsert()
/// (and many other row operations) will produce strange results if used
/// with uninitialized rows. Memory is initialized automatically in the
/// following cases: allocation with ffAlloc(), copying a row with
/// ffCopyRow(), reading a row from a file. A row can be initialized
/// manually by multiplication with zero: ffMulRow(ptr,FF_ZERO,noc).
/// The latter works even if no field has been selected.
/// <br>
///
/// @par Matrices
/// A matrix is stored as a sequence of rows with no extra space between them.
/// The rows and columns of a matrix are numbered 0,1,2,...n-1.
/// Column vectors may be stored as matrices with one column but because of the
/// padding each mark of the vector will occupy a long int.


/// @fn ffAdd(FEL,FEL)
/// Finite field addition.
/// This function returns the sum of two field elements. Before calling ffAdd(), the
/// field must have been selected with ffSetField(). The arguments are not checked. If either
/// argument is not in the allowed range the result is undefined and the program may crash.
/// ffAdd() may be implemented as a macro. In this case, it is guaranteed that both arguments
/// are evaluated exactly once.

/// @fn ffSub(FEL,FEL)
/// Finite field subtraction.
/// This function returns the difference of two field elements. Before calling ffSub(), the
/// field must have been selected with ffSetField(). The arguments are not checked. If either
/// argument is not in the allowed range the result is undefined and the program may crash.
/// ffSub() may be implemented as a macro. In this case, it is guaranteed that both arguments
/// are evaluated exactly once.

/// @fn ffMul(FEL,FEL)
/// Finite field multiplication.
/// This function returns the product of two field elements. Before calling ffMul(), the
/// field must have been selected with ffSetField(). The arguments are not checked. If either
/// argument is not in the allowed range the result is undefined and the program may crash.
/// ffMul() may be implemented as a macro. In this case, it is guaranteed that both arguments
/// are evaluated exactly once.

/// @fn ffDiv(FEL,FEL)
/// Finite field division.
/// This function returns the quotient of two field elements. Before calling ffDiv(), the
/// field must have been selected with ffSetField(). The arguments are not checked. If either
/// argument is not in the allowed range or if the denominator is zero, the result is undefined
/// and the program may crash.
/// ffDiv() may be implemented as a macro. In this case, it is guaranteed that both arguments
/// are evaluated exactly once.

/// @fn ffNeg(FEL)
/// Finite field negative.
/// This function returns the additive inverse a field element. Before calling ffInv(), the
/// field must have been selected with ffSetField(). The argument is not checked. If you pass
/// an invalid value, the result is undefined and the program may crash.
/// ffNeg() may be implemented as a macro. In this case, it is guaranteed that the argument
/// is evaluated exactly once.

/// @fn ffInv(FEL)
/// Finite field inversion.
/// This function returns the multiplicative inverse a field element. Before calling ffInv(), the
/// field must have been selected with ffSetField(). The argument is not checked. If you pass
/// an invalid value or zero, the result is undefined and the program may crash.
/// ffInv() may be implemented as a macro. In this case, it is guaranteed that the argument
/// is evaluated exactly once.

////////////////////////////////////////////////////////////////////////////////////////////////////

static FILE *OpenTableFile(int fl)
{
   char fn[250];
   FILE *fd;

   /* Try to open the table file
      -------------------------- */
   sprintf(fn,"p%3.3d.zzz",fl);
   if ((fd = sysFopen(fn,"rb::lib:noerror")) != NULL) {
      return fd;
   }

   /* Create the table file.
      ---------------------- */
   if (ffMakeTables(fl) != 0) {
      mtxAbort(MTX_HERE,"Unable to build arithmetic tables");
      return NULL;
   }
   fd = sysFopen(fn,"rb::lib");
   return fd;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static int ReadTableFile(FILE *fd, int field)
{
   uint32_t hdr[5];

   /* Read header, perform some checks
      -------------------------------- */
   sysRead32(fd,hdr,5);
   if ((hdr[2] != field) || (hdr[1] < 0) || (hdr[1] > field) ||
       (hdr[0] <= 1) || (hdr[2] % hdr[0] != 0) || (hdr[3] < 1) || (hdr[3] > 8)) {
      mtxAbort(MTX_HERE,"Table file is corrupted");
      return -1;
   }
   ffChar = hdr[0];
   ffGen = (FEL) hdr[1];
   MPB = hdr[3];
   if (hdr[4] != (long) MTX_ZZZVERSION) {
      mtxAbort(MTX_HERE,"Bad table file version: expected %d, found %d",
                 (int)MTX_ZZZVERSION,(int)hdr[4]);
      fclose(fd);
      return -1;
   }

   uint32_t subfields[4];

   /* Read tables
      ----------- */
   sysRead8(fd, mtx_tmult, sizeof(mtx_tmult));
   sysRead8(fd, mtx_tadd, sizeof(mtx_tadd));
   sysRead8(fd, mtx_tffirst,sizeof(mtx_tffirst));
   sysRead8(fd, mtx_textract,sizeof(mtx_textract));
   sysRead8(fd, mtx_taddinv,sizeof(mtx_taddinv));
   sysRead8(fd, mtx_tmultinv,sizeof(mtx_tmultinv));
   sysRead8(fd, mtx_tnull,sizeof(mtx_tnull));
   sysRead8(fd, mtx_tinsert,sizeof(mtx_tinsert));
   sysRead32(fd,subfields,MTX_MAXSUBFIELDS);
   sysRead8(fd,mtx_embed, MTX_MAXSUBFIELDS * MTX_MAXSUBFIELDORD);
   sysRead8(fd,mtx_restrict, MTX_MAXSUBFIELDS * 256);

   // Copy subfields to public table
   memset(mtx_subfields, 0, sizeof(mtx_subfields));
   for (int i = 0; i < 4 && subfields[i] >= 2; ++i) {
      mtx_subfields[i] = (int) subfields[i];
   }

   ffOrder = field;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Sets the field order.
/// This function sets the current field to GF(@em field) and initializes the field arithmetic.
/// Most kernel functions require that a field has been selected before they are used.
/// @param field Field order.
/// @return 0 on success, -1 otherwise.

int ffSetField(int field)
{
   FILE *fd;
   int result;

   if ((field == ffOrder) || (field < 2)) {
      return 0;
   }
   fd = OpenTableFile(field);
   if (fd == NULL) {
      mtxAbort(MTX_HERE,"Cannot open table file for GF(%d)", field);
      return -1;
   }
   const int context = mtxBegin("Loading arithmetic tables for GF(%d)", field);
   result = ReadTableFile(fd,field);
   mtxEnd(context);
   fclose(fd);
   return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Long ints per row.
static int lpr(int noc)
{
   const int mpl = MPB * sizeof(long);
   return (noc + mpl - 1) / mpl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Calculate row size.
/// Returns the number of bytes occupied in memory by a row of @em noc Elements.
/// The row size is always a multiple of <tt>sizeof(long)</tt>. Depending on the number of
/// columns there may be unused padding bytes at the end of the row.

size_t ffRowSize(int noc)
{
   MTX_ASSERT(noc >= 0);
   return lpr(noc) * sizeof(long);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the number of bytes occupied in memory by @a nor rows of @a noc elements.
/// @a nor may be negative, resulting in a negative return value such that
/// <tt>ffSize(nor,-noc) == -ffSize(nor,noc)</tt>. Thus, @c ffSize() can be used to calculate
/// row pointer differences in both directions.
/// @a noc must be greater than or equal to zero.

ssize_t ffSize(int nor, int noc)
{
   return nor == 0 ? 0 : nor * ffRowSize(noc);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Number of used bytes in a row.
/// This function returns the number of bytes that are actually used by a row of @em noc Elements,
/// i.e., without counting the padding bytes. This number is less than or equal to
/// <tt>ffRowSize(noc)</tt>.

size_t ffRowSizeUsed(int noc)
{
   if (noc == 0) {
      return 0;
   }
   MTX_ASSERT(noc > 0);
   return (noc - 1) / MPB + 1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Embed an element of a subfield.
/// @param a Element of the subfield.
/// @param subfield Subfield order. Must be a divisor of the current field order.
/// @return @em a, embedded into the current field.

FEL ffEmbed(FEL a, int subfield)
{
   int i;

   if (subfield == ffOrder) {
      return a;
   }
   for (i = 0; i < 4; ++i) {
      if (mtx_subfields[i] == subfield) {
         if (a >= subfield) {
	    mtxAbort(MTX_HERE,"Invalid field element %d in GF(%d)",(int) a, subfield);
	 }
         return mtx_embed[i][a];
      }
   }
   mtxAbort(MTX_HERE,"Cannot embed GF(%d) into GF(%d)",(int)subfield,(int)ffOrder);
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Restrict a field element to a subfield. The returned value can be used after switching to
/// the subfield.
///
/// The function fails (aborting the program) if
/// * \a subfield is not an integer root of the current field order, or
/// * the element \a a is not a member of the subfield of order \a subfield.

FEL ffRestrict(FEL a, int subfield)
{
   int i;

   if (subfield == ffOrder) {
      return a;
   }
   for (i = 0; i < 4; ++i) {
      if (mtx_subfields[i] == subfield) {
         FEL result = mtx_restrict[i][a];
         if (result > subfield)
            mtxAbort(MTX_HERE, "Field element is not in GF(%d) < GF(%d)", subfield, ffOrder);
         return result;
      }
   }
   mtxAbort(MTX_HERE,"Cannot restrict from GF(%d) to GF(%d)",ffOrder, subfield);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Add two rows.
/// This function adds src to dest. Field order and row size must have been set before.
/// @param dest The row to add to.
/// @param src The row to add.
/// @return Always returns dest.

PTR ffAddRow(PTR dest, PTR src, int noc)
{
   register int i;

   if (ffChar == 2) {   /* characteristic 2 is simple... */
      register long *l1 = (long *) dest;
      register long *l2 = (long *) src;
      for (i = lpr(noc); i != 0; --i) {
         register long x = *l2++;
         if (x != 0) { *l1 ^= x; }
         l1++;
      }
   } else {             /* any other characteristic */
      register BYTE *p1 = dest;
      register BYTE *p2 = src;
      for (i = ffRowSize(noc); i != 0; --i) {
         register int x = *p2++;
         if (x != 0) { *p1 = mtx_tadd[*p1][x]; }
         p1++;
      }
   }
   return dest;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Adds a row to another row, starting at the given column.
/// This is an optimized version of @ref ffAddRow to be used in row cleaning operations.
/// The function assumes that both source and target rows have already been partially cleaned
/// and contain only zeroes before @a firstcol. If this is not the case, the result is
/// unspecified.
///
/// @param dest The row to add to.
/// @param src The row to add.
/// @param first First column to add.
/// @param noc Row size (number of columns).

void ffAddRowPartial(PTR dest, PTR src, int first, int noc)
{
   MTX_ASSERT(first < noc);

   if (ffChar == 2)     /* characteristic 2 is simple... */
   {
      const int firstl = first / MPB / sizeof(long);
      long *l1 = (long *) dest + firstl;
      long *l2 = (long *) src + firstl;
      for (int i = ffRowSize(noc) / sizeof(long) - firstl; i > 0; --i) {
         long x = *l2++;
         *l1 ^= x;
         l1++;
      }
   }
   else {               /* any other characteristic */
      BYTE *p1 = dest + first / MPB;
      BYTE *p2 = src + first / MPB;
      for (int i = ffRowSize(noc) - first / MPB; i != 0; --i) {
         int x = *p2++;
         *p1 = mtx_tadd[*p1][x];
         p1++;
      }
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Multiply a vector by a field element.
/// This function multiplies each element of @a row by @a mark.
/// The field order must have been set before.
///
/// Multiplying a row with zero (FF_ZERO) initializes all elements to zero and is permitted even if
/// @a row points to uninitialized memory. Furthermore, multiplying with FF_ZERO fills unused bytes
/// at the end of the row with zeroes (which may be different from FF_ERO.

void ffMulRow(PTR row, FEL mark, int noc)
{
   MTX_ASSERT_DEBUG(isFel(mark));
   if (mark == FF_ZERO) {
      memset(row, 0, ffRowSize(noc));
   } else if (mark != FF_ONE) {
      const uint8_t* const multab = mtx_tmult[mark];
      uint8_t* m = row;
      for (int i = ffRowSize(noc); i != 0; --i) {
         if (*m != 0) 
            *m = multab[*m];
         ++m;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Adds a multiple of a row (@a src) to another row (@a dest). Both rows must have the same
/// size (@a noc).

void ffAddMulRow(PTR dest, PTR src, FEL f, int noc)
{
   MTX_ASSERT_DEBUG(isFel(f));
   if (f == FF_ONE) {
      ffAddRow(dest,src, noc);
   }
   else if (f != FF_ZERO) {
      uint8_t *multab = mtx_tmult[f];
      uint8_t *p1 = dest;
      uint8_t *p2 = src;
      for (int i = ffRowSize(noc); i != 0; --i) {
         if (*p2 != 0) {
            *p1 = mtx_tadd[*p1][multab[*p2]];
         }
         ++p1;
         ++p2;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Adds a multiple of a row, starting at the given column.
///
/// This is an optimized version of @ref ffAddMulRow to be used in row cleaning operations.
/// The function assumes that both source and target rows have already been partially cleaned
/// and contain only zeroes before @a firstcol. If this is not the case, the result is
/// unspecified.
///
/// @param dest The row to add to.
/// @param src The row to add.
/// @param f The multiplier.
/// @param firstcol First column to add.
/// @param noc Row size (number of columns).

void ffAddMulRowPartial(PTR dest, PTR src, FEL f, int firstcol, int noc)
{
   MTX_ASSERT_DEBUG(isFel(f));
   MTX_ASSERT_DEBUG(firstcol >= 0 && firstcol < noc);

   if (f == FF_ONE) {
      ffAddRowPartial(dest, src, firstcol, noc);
   }
   else if (f != FF_ZERO) {
      BYTE * const multab = mtx_tmult[f];
      BYTE *p1 = dest + firstcol / MPB;
      BYTE *p2 = src + firstcol / MPB;
      for (int i = ffRowSize(noc) - firstcol / MPB; i != 0; --i) {
         if (*p2 != 0) {
            *p1 = mtx_tadd[*p1][multab[*p2]];
         }
         ++p1;
         ++p2;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Convert integer to field element.
/// This function, together with ffFromInt(), defines a bijection between field elements and
/// the set of integers {0, 1, ... q-1}, where q is the field order. This mapping is used when
/// field elements are converted to or from text form and has the following properties:
/// - ffFromInt(0) if the zero element
/// - ffFromInt(1) if the unit element
/// - The restriction to {0,..,p-1} (with the usual arithmetic mod p) is an isomorphism of
///   Z/(pZ) to the prime field.

FEL ffFromInt(int l)
{
   register int f;
   f = l % ffOrder;
   if (f < 0) { f += ffOrder; }
   return (FEL) f;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int ffToInt(FEL f)
{
   return (int) f;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Multiply a vector by a matrix.
/// This function multiplies the vector @a row from the right by the matrix @a mat and
/// stores the result into @a result.
///
/// @attention @em result and @em row must not overlap. Otherwise the result is undefined.
///
/// @param row The source vector (\a nor columns).
/// @param matrix The matrix (\a nor by \a noc).
/// @param nor number of rows in the matrix.
/// @param noc number of columns in the matrix.
/// @param[out] result The resulting vector (@a noc columns).

void ffMapRow(PTR row, PTR matrix, int nor, int noc, PTR result)
{
   register int i;
   register FEL f;
   BYTE *m = (BYTE *) matrix;

   // Check that result and row do not overlap.
   MTX_ASSERT(ffGetPtr(row, 1, nor) <= result || row >= ffGetPtr(result, 1, noc));

   // Initialize result.
   ffMulRow(result, FF_ZERO, noc);

   const int LPR = lpr(noc);
   if (ffOrder == 2) {      
      // Use bitwise operations
      register long *x1 = (long *) matrix;
      register BYTE *r = (BYTE *) row;

      for (i = nor; i > 0; ++r) {
         register BYTE mask;
         if (*r == 0) {
            i -= 8;
            x1 += 8 * LPR;
            continue;
         }
         for (mask = 0x80; mask != 0 && i > 0; mask >>= 1, --i) {
            if ((mask & *r) == 0) {
               x1 += LPR;       // Skip this row
            }
            else
            {
               register long *x2 = (long *)result;
               register int k;
               for (k = LPR; k; --k) {
                  *x2++ ^= *x1++;
               }
            }
         }
      }
   } else {             /* Any other field */
      const size_t rowSize = ffRowSize(noc);
      register BYTE *brow = (BYTE *) row;
      register int pos = 0;
      for (i = nor; i > 0; --i) {
         f = mtx_textract[pos][*brow];
         if (++pos == (int) MPB) {
            pos = 0;
            ++brow;
         }
         if (f != FF_ZERO) {
            register BYTE *v = m;
            register BYTE *r = result;
            register int k = rowSize;
            if (f == FF_ONE) {
               for (; k != 0; --k) {
                  if (*v != 0) {
                     *r = mtx_tadd[*r][*v];
                  }
                  ++r;
                  ++v;
               }
            } else {
               register BYTE *multab = mtx_tmult[f];
               for (; k != 0; --k) {
                  if (*v != 0) {
                     *r = mtx_tadd[multab[*v]][*r];
                  }
                  ++v;
                  ++r;
               }
            }
         }
         m += rowSize;                 /* next row */
      }
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Scalar Product of Two Vectors.
/// Given two vectors @f$a=(a_i)@f$ and @f$b=(b_i)@f$, this function calculates the
/// scalar product @f$p=\sum_ia_ib_i@f$.
/// @param a The first vector.
/// @param b The second vector.
/// @return Scalar product of the two vectors.

FEL ffScalarProduct(PTR a, PTR b, int noc)
{
   const uint8_t *ap = (const uint8_t *) a;
   const uint8_t *bp = (const uint8_t *) b;
   FEL f = FF_ZERO;

   int i = noc;
   // full bytes
   for (; i >= MPB; i -= MPB) {
      if (*ap != 0 && *bp != 0) {
         for (int k = 0; k < MPB; ++k) {
            f = ffAdd(f,ffMul(mtx_textract[k][*ap],mtx_textract[k][*bp]));
         }
      }
      ++ap;
      ++bp;
   }
   // last byte may have less than MBP marks
   while (i > 0) {
      --i;
      f = ffAdd(f,ffMul(mtx_textract[i][*ap],mtx_textract[i][*bp]));
   }

   return f;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Extract one column of a matrix.
/// This function extracts one column out of a matrix and converts it into a row vector.
/// in the call. The result is a row with @a nor entries, and the output buffer,
/// @a result, must be large enough to store a vector of this size.
/// @a mat and @a result must not overlap. If they do, the result is undefined.
///
/// @param mat Pointer to the matrix (@a nor by @a noc).
/// @param nor Number of rows in matrix.
/// @param nor Number of columns in matrix.
/// @param col Column to extract (starting with 0).
/// @param result Pointer to buffer for the extracted column (@a nor columns).

void ffExtractColumn(PTR mat, int nor, int noc, int col, PTR result)
{
   register BYTE *x = (BYTE *)mat + (col / MPB);
   register BYTE *extab = mtx_textract[col % MPB];
   register BYTE a = 0;
   register int ind = 0;
   register BYTE *y = result;
   register int count;

   for (count = nor; count > 0; --count) {
      a = (BYTE) (a + mtx_tinsert[ind][extab[*x]]);
      if (++ind == MPB) {
         *y++ = a;
         a = 0;
         ind = 0;
      }
      x += ffRowSize(noc);
   }
   if (ind != 0) { *y = a; }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Insert a mark into a row
/// This function inserts the field element @em mark at position @em col into @em row.
/// Column indexes start with 0.
/// Before this function can be used, the field must be selected with ffSetField().
///
/// @note @c ffInsert() fails on negative column numbvers but does not detect writing beyond
/// the end of @a row. Doing so results in undefined behaviour.
///
/// @param row Pointer to the row.
/// @param col Insert position (0-based).
/// @param mark Value to insert.

void ffInsert(PTR row, int col, FEL mark)
{
   MTX_ASSERT(col >= 0);
   MTX_ASSERT_DEBUG(isFel(mark));

   register BYTE *loc = (BYTE *)row + (col / MPB);
   register int idx = col % MPB;
   *loc = (BYTE) (mtx_tnull[idx][*loc] + mtx_tinsert[idx][mark]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Extract a mark from a row
/// This function returns the entry at position @a col of @a row.
/// Note that column indexes start with 0, i.e., ffExtract(row,0) returns the first entry.
/// Like ffInsert(), this function does not depend on the current row size.
/// The index, @a col, is not checked. Reading at negative positions or beyond the end of the row
/// results in undefined behaviour.
///
/// @param row Pointer to the row.
/// @param col Index of mark to extract (0 based).
/// @return The specified mark
///
/// @see FfInsert

FEL ffExtract(PTR row, int col)
{
   FEL result = mtx_textract[col % MPB][((BYTE *)row)[col / MPB]];
   MTX_ASSERT_DEBUG(isFel(result));
   return result;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Find pivot column.
/// This function finds the first non-zero mark in a row vector.
/// The mark is stored into <tt>*mark</tt> and its position (counting from 0) is returned.
/// If the whole vector is zero, %ffFindPivot()
/// returns -1 and leaves <tt>*mark</tt> unchanged.
/// @param row Pointer to the row vector (@a noc columns).
/// @param mark Buffer for pivot element.
/// @param noc Number of columns.
/// @return Index of the first non-zero entry in @a row or -1 if all entries are zero.

int ffFindPivot(PTR row, FEL *mark, int noc)
{
   register long *l = (long *) row;
   register int idx;
   register BYTE *m;

   const int LPR = lpr(noc);
   for (idx = 0; idx < LPR && *l == 0; ++idx, ++l) {
   }
   if (idx == LPR) {
      return -1;        // all zero
   }
   idx = idx * sizeof(long) * MPB;
   m = (BYTE *)l;
   while (*m == 0) {
      ++m;
      idx += MPB;
   }
   idx += mtx_tffirst[*m][1];
   if (idx >= noc) {          // Ignore garbage in padding bytes
      return -1;
   }
   *mark = mtx_tffirst[*m][0];
   return idx;
}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
