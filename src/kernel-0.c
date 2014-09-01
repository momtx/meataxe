/* ============================= C MeatAxe ==================================
   File:        $Id: kernel-0.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Finite field arithmetic and common functions. `Small' version
		for field orders q <= 256. Originally based on the `hprout.c'
		written by Klaus Lux.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   Modified by Tom Hoffman, Mathematics Department, University of
   Arizona, USA <hoffmant@math.arizona.edu>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include <string.h>
#include <stdlib.h>
#include "meataxe.h"


/* ------------------------------------------------------------------
   Gobal data
   ------------------------------------------------------------------ */

MTX_DEFINE_FILE_INFO

typedef unsigned char BYTE;
static int MPB = 0;		/* No. of marks per byte */
static int LPR = 0;		/* Long ints per row */





/* ------------------------------------------------------------------
   Argument checking (for debbugging)
   ------------------------------------------------------------------ */

#if defined(_DEBUG)
#define CHECKRANGE(x,lo,hi) if ((x)<(lo)||(x)>(hi)) {\
	MtxError(&Mtx_ThisFile,__LINE__,\
	"ZZZ RANGE CHECK ERROR: %d <= %d <= %d ?\n",\
	(int)(lo),(int)(x),(int)(hi));}
#else
#define CHECKRANGE(x,lo,hi)
#endif

#define CHECKFEL(x) CHECKRANGE(x,0,FfOrder-1)
#define CHECKCOL(x)  CHECKRANGE(x,1,FfNoc)

/**
 ** @defgroup ff Finite Fields
 ** 
 ** The finite field part of the kernel provides finite field arithmetic and
 ** basic operations with vectors and matrices over finite fields.
 ** The kernel cannot operate simultaneously with different finite fields,
 ** because there is a global row size and a global field order which must
 ** be maintained by the higher layers.
 **
 ** There are two finite field modules available: one for small fields (up
 ** to 256) and one for larger fields (up to 2<sup>16</sup>). The finite field
 ** module is selected at compile time.
 ** 
 ** @section ff_dt Basic data types
 ** The kernel defines two basic data types:
 ** @par FEL
 ** represents a single field element
 ** @par PTR
 ** is a pointer to a row vector. PTR may be defined as <tt>*FEL</tt>, but this is
 ** not mandatory.
 **
 ** The kernel also defines two constants:
 ** FF_ZERO is the zero element of the current field, and
 ** FF_ONE is unit element of the current field.
 ** Depending on which kernel you are using, FF_ZERO and FF_ONE need
 ** not be constants.
 ** They may be defined as variables or even function calls.
 ** 
 ** @section ff_intrep Internal data representation
 ** @par 'Small' Kernel (q≤256)
 ** In the "small" kernel, field elements of GF(q) are represented
 ** by the numbers 0,1,...,q-1. The field is defined by its
 ** Conway polynomial p(x), a polynomial of degree n over ℤ<sub>p</sub>[x],
 ** where q=p<sup>n</sup>. Thus, we have a one-to-one correspondence of field
 ** elements a∈GF(q) and polynomials f<sub>a</sub>(x)∈ ℤ<sub>p</sub>[x] of
 
 ** degree less than n. By treating ℤ<sub>p</sub> as a subset of ℤ ---
 ** actually, on the computer, elements of ℤ<sub>p</sub> are represented
 ** by integers --- this is also a polynomial over *ℤ. Now,
 ** calculate f_a(p) giving the number of the field element a.
 ** It follows that the elements of the prime field are represented by
 ** 0,1,...p-1. The number 0 represents the zero element,
 ** and 1 represents the unit element.
 **
 ** @par 'Big' Kernel (q≤65536)
 ** The big version stores field elements in 16-bit integers, i.e., each field element
 ** occupies two bytes. Non-zero elements are stored as their logarithms with
 ** respect to a fixed generator. In particular, the unit element is represented by
 ** the integer 0. The zero element is represented by the special value 0xFFFF.
 **
 ** As a consequence of the different representations of field elements
 ** in the small and big version, there are some rules which should be
 ** respected by all programs:
 ** - Never assign numbers to variables of type FEL or pass
 **   numbers to functions expecting an argument of type FEL.
 ** - Never perform integer arithmetic on variables of type FEL.
 ** - Never printf() or scanf() variables of type FEL.
 ** - Never use the literals 0 or 1 where the zero and unit element
 **   of the field is intended. Instead you should use the FF_ZERO
 **   and FF_ONE constants define in |meataxe.h|.
 ** - Do not cast an integer to FEL or vice versa. Use
 **   FfFromInt() and FfToInt() instead.
 ** 
 ** @section ff_convert Converting between finite field elements and integers
 ** The MeatAxe defines a standard numbering of field elements,
 ** i.e., a bijection between GF(q) and the set {0,1,..,q-1}.
 ** For prime fields, the mapping is defined by assigning the integer 1
 ** to the unit element of the field.
 ** For non-prime fields or forder q=p<sup>n</sup>, each a∈GF(q) is represented
 ** -- modulo the ideal generated by the Conway polynomial of degree n -- by a
 ** unique polynomial f<sub>a</sub>(x)∈GF(p)[x] with deg(f<sub>a</sub>)<n.
 ** Using the standard embedding of GF(p) into ℤ
 ** we can consider f<sub>a</sub> as a polynomial over ℤ.
 ** Then, the number assigned to a is f<sub>a</sub>(p).
 **
 ** The mapping between GF(q) and {0,1,...,q} is provided by two
 ** functions, FfFromInt() and FfToInt(). Since the actual numeric
 ** representation of field elements depends on the kernel, you cannot
 ** convert a FEL to an integer simply by type casting.
 ** 
 ** @section ff_embed Embedding of subfields
 ** In the MeatAxe there is a `standard' generator for each finite
 ** field. The generator for the field currently in use is available in
 ** the FfGen variable. Thus, if a and a' are the MeatAxe
 ** generators of GF(q) and GF(q'), respectively, and q'=q<sup>n</sup>, there
 ** is a standard embedding of GF(q) into GF(q') defined by
 ** $a↪a'<sup>n</sup>. However, field elements which are identified
 ** under this embedding are usually not represented by the same number.
 ** For this reason there are two functions, FfEmbed() and
 ** FfRestrict(), which provide the embedding of subfields into the
 ** current field. Note that the MeatAxe is not well suited for
 ** calculations involving different fields at the same time because the
 ** arithmetic uses lookup tables which are specific to each field
 **
 ** Here is a short example. The following code converts a vector over GF(3)
 ** to GF(27). The MeatAxe cannot handle two fields at the same time, so
 ** it is necessary to unpack the row over GF(3), change to GF(27) and pack
 ** the embedded elements into a new row.
 ** @code
 ** PTR row1, row2;
 ** FEL buf[10];
 ** ...
 ** FfSetField(10);
 ** for (i = 0; i &lt; 10; ++i) buf[i] = FfExtract(row1,i);
 ** FfSetField(27);
 ** for (i = 0; i &lt; 10; ++i) FfInsert(row2,i,FfEmbed(buf[i],3));
 ** @endcode
 **/

/**
 ** @defgroup ff2 Other Kernel Functions
 **/

/**
 ** @defgroup ffrow Row Operations
 ** The fuctions below perform operations on row vectors.
 ** All row operations at kernel level use a global row size which can be
 ** set with FfSetNoc(). The current row size is available in the global
 ** variable FfNoc.
 **
 ** @par Internal data format
 **
 ** A (row) vector is stored in memory as an array of bytes. The actual
 ** size of a row depends on the number of marks in the row and on the
 ** field order, because several field elements may be packed into one
 ** byte. However, the size is always a multiple of sizeof(long).
 ** Thus, there may be unused bytes at the end of a row. The contents of
 ** these extra bytes (and even the contents of unused bits in partially
 ** used bytes) is not defined. For this reason, memory must be initialized
 ** before it is used, or else some functions as  FfCmpRows() may fail.
 **
 ** Packing of field elements is used for field orders less than 17. Let q be
 ** the field order and m the largest natural number with q<sup>m</sup>≤256.
 ** Then, m elements k<sub>0</sub>,...k<sub>n-1</sub> are packed into one byte
 ** as k<sub>0</sub>+k<sub>1</sub>q+k<sub>2</sub>q<sup>2</sup>+...
 ** Packing of field elements is used exclusively for rows (vectors),
 ** not for polynomials or other data types.
 **
 ** Because of packing you cannot access the marks of a row with the usual bracket notation.
 ** Instead, the MeatAxe library provides functions which store and extract marks.
 **
 ** Packed rows must be initialized before they are used. FfInsert()
 ** (and many other row operations) will produce strange results if used
 ** with uninitialized rows. Memory is initialized automatically in the
 ** following cases: allocation with FfAlloc(), copying a row with
 ** FfCopyRow(), reading a row from a file. A row can be initialized
 ** manually by multiplication with zero: FfMulRow(ptr,FF_ZERO).
 ** The latter works even if no field has been selected.
 **
 ** A matrix is stored as a sequence of rows with no extra space between them.
 ** The rows and columns of a matrix are numbered 0,1,2,...n-1.
 ** Column vectors may be stored as matrices with one column but each mark
 ** of the vector will occupy 4 Bytes.
 **/


/**
 ** @addtogroup ff
 ** @{
 **/

/**
 ** Row size.
 ** This variable contains the size of a single row in memory. Its value is always equal to
 ** <tt>FfRowSize(FfNoc)</tt>. The row size os always a multiple of sizeof(long).
 **/

size_t FfCurrentRowSize = (size_t) -1;

/**
 ** @ingroup ff
 ** I/O row size.
 ** This variable contains the number of bytes occupied by a row when stored in a data
 ** file. Its value is always equal to <tt>FfTrueRowSize(FfNoc)</tt>. Since there is no
 ** padding in data files, FfCurrentRowSizeIo is usually smaller than FfCurrentRowSize.
 **/

int FfCurrentRowSizeIo = -1;


/**
 ** @fn FfAdd(FEL,FEL)
 ** Finite field addition.
 ** This function returns the sum of two field elements. Before calling FfAdd(), the
 ** field must have been selected with FfSetField(). The arguments are not checked. If either
 ** argument is not in the allowed range the result is undefined and the program may crash.
 ** FfAdd() may be implemented as a macro. In this case, it is guaranteed that both arguments
 ** are evaluated exactly once.
 **/

/**
 ** @fn FfSub(FEL,FEL)
 ** Finite field subtraction.
 ** This function returns the difference of two field elements. Before calling FfSub(), the
 ** field must have been selected with FfSetField(). The arguments are not checked. If either
 ** argument is not in the allowed range the result is undefined and the program may crash.
 ** FfSub() may be implemented as a macro. In this case, it is guaranteed that both arguments
 ** are evaluated exactly once.
 **/

/**
 ** @fn FfMul(FEL,FEL)
 ** Finite field multiplication.
 ** This function returns the product of two field elements. Before calling FfMul(), the
 ** field must have been selected with FfSetField(). The arguments are not checked. If either
 ** argument is not in the allowed range the result is undefined and the program may crash.
 ** FfMul() may be implemented as a macro. In this case, it is guaranteed that both arguments
 ** are evaluated exactly once.
 **/

/**
 ** @fn FfDiv(FEL,FEL)
 ** Finite field division.
 ** This function returns the quotient of two field elements. Before calling FfDiv(), the
 ** field must have been selected with FfSetField(). The arguments are not checked. If either
 ** argument is not in the allowed range or if the denominator is zero, the result is undefined
 ** and the program may crash.
 ** FfDiv() may be implemented as a macro. In this case, it is guaranteed that both arguments
 ** are evaluated exactly once.
 **/

/**
 ** @fn FfNeg(FEL)
 ** Finite field negative.
 ** This function returns the additive inverse a field element. Before calling FfInv(), the
 ** field must have been selected with FfSetField(). The argument is not checked. If you pass
 ** an invalid value, the result is undefined and the program may crash.
 ** FfNeg() may be implemented as a macro. In this case, it is guaranteed that the argument
 ** is evaluated exactly once.
 **/

/**
 ** @fn FfInv(FEL)
 ** Finite field inversion.
 ** This function returns the multiplicative inverse a field element. Before calling FfInv(), the
 ** field must have been selected with FfSetField(). The argument is not checked. If you pass
 ** an invalid value or zero, the result is undefined and the program may crash.
 ** FfInv() may be implemented as a macro. In this case, it is guaranteed that the argument
 ** is evaluated exactly once.
 **/


/* -----------------------------------------------------------------
   OpenTableFile()
   ------------------------------------------------------------------ */

static FILE *OpenTableFile(int fl)
{
    char fn[250];
    FILE *fd;

    /* Try to open the table file
       -------------------------- */
    sprintf(fn,"p%3.3d.zzz",fl);
    if ((fd = SysFopen(fn,FM_READ|FM_LIB|FM_NOERROR)) != NULL)
	return fd;

    /* Create the table file.
       ---------------------- */
    if (FfMakeTables(fl) != 0)
	MTX_ERROR("Unable to build arithmetic tables");
    fd = SysFopen(fn,FM_READ|FM_LIB);
    return fd;
}




/* -------------------------------------------------------------------------
   ReadTableFile()
   -------------------------------------------------------------------------- */

static int ReadTableFile(FILE *fd, int field)
{
    long hdr[5];

    /* Read header, perform some checks
       -------------------------------- */
    if (SysReadLong(fd,hdr,5) != 5)
    {
	MTX_ERROR("Cannot read table file header");
	return -1;
    }
    if (hdr[2] != field || hdr[1] < 0 || hdr[1] > field ||
        hdr[0] <= 1 || hdr[2] % hdr[0] != 0 || hdr[3] < 1 || hdr[3] > 8)
    {	
	MTX_ERROR("Table file is corrupted");
	return -1;
    }
    FfChar = hdr[0];
    FfGen = (FEL) hdr[1];
    MPB = hdr[3];
    if (hdr[4] != (long) ZZZVERSION)
    {
	MTX_ERROR2("Bad table file version: expected %d, found %d",
	    (int)ZZZVERSION,(int)hdr[4]);
	fclose(fd);
	return -1;
    }

    /* Read tables
       ----------- */
    if (fread(mtx_tmult,4,sizeof(mtx_tmult)/4,fd) != sizeof(mtx_tmult)/4 ||
        fread(mtx_tadd,4,sizeof(mtx_tadd)/4,fd) != sizeof(mtx_tadd)/4 ||
        fread(mtx_tffirst,1,sizeof(mtx_tffirst),fd) != sizeof(mtx_tffirst) ||
        fread(mtx_textract,1,sizeof(mtx_textract),fd) != sizeof(mtx_textract) ||
        fread(mtx_taddinv,1,sizeof(mtx_taddinv),fd) != sizeof(mtx_taddinv) ||
        fread(mtx_tmultinv,1,sizeof(mtx_tmultinv),fd) != sizeof(mtx_tmultinv) ||
        fread(mtx_tnull,1,sizeof(mtx_tnull),fd) != sizeof(mtx_tnull) ||
        fread(mtx_tinsert,1,sizeof(mtx_tinsert),fd) != sizeof(mtx_tinsert) ||
        SysReadLong(fd,mtx_embedord,4) != 4 ||
        fread(mtx_embed,16,4,fd) != 4 ||
        fread(mtx_restrict,256,4,fd) != 4
       )
    {
        MTX_ERROR("Error reading table file");
	return -1;
    }
    FfOrder = field;
    FfSetNoc(FfOrder);
    return 0;
}






/**
 ** Set the field order.
 ** This function sets the current field to GF(@em field) and initializes the field arithmetic.
 ** Most kernel functions require that a field has been selected before they are used.
 ** @param field Field order.
 ** @return 0 on success, -1 otherwise.
 **/

int FfSetField(int field)

{
    FILE *fd;
    int result;

    if (field == FfOrder || field < 2)
	return 0;
    fd = OpenTableFile(field);
    if (fd == NULL)
    {
	MTX_ERROR1("Cannot open table file for GF(%d)",(int)field);
	return -1;
    }
    result = ReadTableFile(fd,field);
    fclose(fd);
    return result;
}


/**
 ** Set the row length.
 ** This function sets the current row size, which is used for low-level row operations
 ** such as FfAddRow().
 ** @param noc Number of columns.
 ** @return 0 on success, -1 otherwise.
 **/

int FfSetNoc(int noc)
{
    MTX_ASSERT(noc >= 0);
    FfNoc = noc;
    if (noc == 0)
    {
	LPR = 0;
	FfCurrentRowSize = 0;
	FfCurrentRowSizeIo = 0;
    }
    else
    {
	LPR = (noc-1) / (MPB * sizeof(long)) + 1;	/* long's per row */
#ifdef ASM_MMX
	if (LPR % 2 != 0)
	    ++LPR;
#endif
	FfCurrentRowSize = LPR * sizeof(long);
        FfCurrentRowSizeIo = (noc - 1) / MPB + 1;
    }
    return 0;
}


/**
 ** Calculate row size.
 ** Returns the number of bytes occupied in memory by a row of @em noc Elements.
 ** The row size is always a multiple of <tt>sizeof(long)</tt>. Depending on the number of
 ** columns there may be unused padding bytes at the end of the row.
 **/

size_t FfRowSize(int noc)
{
    if (noc == 0)
	return 0;
    MTX_ASSERT(noc > 0);
    return ((noc-1) / (MPB * sizeof(long)) + 1) * sizeof(long);
}


/**
 ** Number of used bytes in a row.
 ** This function returns the number of bytes that are actually used by a row of @em noc Elements,
 ** i.e., without counting the padding bytes. This number is less than or equal to
 ** <tt>FfRowSize(noc)</tt>.
 **/

size_t FfTrueRowSize(int noc)

{
    if (noc == 0)
	return 0;
    MTX_ASSERT(noc > 0);
    return (noc-1) / MPB + 1;
}




/**
 ** Embed a subfield.
 ** @param a Element of the subfield field.
 ** @param subfield Subfield order. Must be a divisor of the current field order.
 ** @return @em a, embedded into the current field.
 **/ 

FEL FfEmbed(FEL a, int subfield)

{   int i;

    if (subfield == FfOrder)
	return a;
    for (i = 0; mtx_embedord[i] != subfield && i < 4; ++i);
    if (i >= 4)
	MTX_ERROR2("Cannot embed GF(%d) into GF(%d)",(int)subfield,(int)FfOrder);
    return mtx_embed[i][a];
}



/**
 ** Restrict to a subfield.
 ** This function restricts a field element from the current field to a subfield. The return
 ** value represents the same element as @em a but with respect to the subfield. In general,
 ** the element has a different integer representation in the subfield. Consequently, you cannot
 ** use the return value for field arithmetic until you change to the subfield with
 ** Of course, the argument must be an element of the subfield. Otherwise the result is undefined.
 ** <tt>FfSetField(subfield)</tt>.
 ** @param a Element of the current field.
 ** @param subfield Subfield order. Must be a divisor of the current field order.
 **/

FEL FfRestrict(FEL a, int subfield)
{
    int i;

    if (subfield == FfOrder)
	return a;
    for (i = 0; mtx_embedord[i] != subfield && i < 4; ++i);
    if (i >= 4)
    {
	MTX_ERROR2("Cannot restrict GF(%d) to GF(%d)",(int)FfOrder,
	    (int)subfield);
    }
    return mtx_restrict[i][a];
}



#ifdef ASM_MMX
static void __inline__ FastXor(void *dest, const void *src)
{
__asm__(
	
	"    pushl %ebx\n"
	"    pushl %ecx\n"
	"    pushl %edx\n"

	"    movl 8(%ebp),%ecx\n"
        "    movl 12(%ebp),%ebx\n"
        "    movl LPR,%edx\n"
        "    sarl $1,%edx\n"
        "    je .FASTXOR_1\n"
        "    .align 16\n"
	".FASTXOR_2:\n"
        "    movq (%ebx),%mm0\n"
        "    addl $8,%ebx\n"
        "    pxor (%ecx),%mm0\n"
        "    movq %mm0,(%ecx)\n"
        "    addl $8,%ecx\n"
        "    decl %edx\n"
        "    jne .FASTXOR_2\n"
	".FASTXOR_1:\n"
	"    popl %edx\n"
	"    popl %ecx\n"
	"    popl %ebx\n"
	);
}
#endif


/**
 ** @}
 **/


/**
 ** @addtogroup ffrow
 ** @{
 **/

/**
 ** Add two rows.
 ** This function adds src to dest. Field order and row size must have been set before.
 ** @param dest The row to add to.
 ** @param src The row to add.
 ** @return Always returns dest.
 **/

PTR FfAddRow(PTR dest, PTR src)
{
    register int i;

    if (FfChar == 2)	/* characteristic 2 is simple... */
    {	
#ifdef ASM_MMX
    /* This assumes Intel with 4 bytes per long, but MMX implies Intel anyway.*/
	__asm__(
	"    pushl %ebx\n"
	"    pushl %ecx\n"
	"    pushl %edx\n"

	"    movl 8(%ebp),%ecx\n"
        "    movl 12(%ebp),%ebx\n"
        "    movl LPR,%edx\n"
        "    sarl $1,%edx\n"
        "    je .ADDROW2\n"
        "    .align 16\n"
	".ADDROW1:\n"
        "    movq (%ebx),%mm0\n"
        "    addl $8,%ebx\n"
        "    pxor (%ecx),%mm0\n"
        "    movq %mm0,(%ecx)\n"
        "    addl $8,%ecx\n"
        "    decl %edx\n"
        "    jne .ADDROW1\n"
	".ADDROW2:\n"
	"    popl %edx\n"
	"    popl %ecx\n"
	"    popl %ebx\n"
	);
#else
	register long *l1 = (long *) dest;
	register long *l2 = (long *) src;
	for (i = LPR; i != 0; --i)
	{
	    register long x = *l2++;
	    if (x != 0) *l1 ^= x;
	    l1++;
	}
#endif
    }
    else		/* any other characteristic */
    {
#ifdef ASM_MMX
        register BYTE *p1 = dest;
        register unsigned long *p2 = (unsigned long *) src;
        for (i = LPR; i != 0; --i)
        {
            register unsigned long a;
            if ((a = *p2++) != 0) {
                *p1++ = mtx_tadd[*p1][a & 0xffL];
                a >>= 8;
                *p1++ = mtx_tadd[*p1][a & 0xffL];
                a >>= 8;
                *p1++ = mtx_tadd[*p1][a & 0xffL];
                a >>= 8;
                *p1++ = mtx_tadd[*p1][a & 0xffL];
            } else
              p1 += 4;
        }
#else
	register BYTE *p1 = dest;
	register BYTE *p2 = src;
	for (i = FfTrueRowSize(FfNoc); i != 0; --i)
	{
	    register int x = *p2++;
	    if (x != 0) *p1 = mtx_tadd[*p1][x];
	    p1++;
	}
#endif
    }
    return dest;
}



/**
 ** Add a part two rows.
 ** This works like FfAddRow(), but the operation is performed only on a given range of
 ** columns. Note that the working range is not specified as column indexes but in units of
 ** long integers!
 ** @param dest The row to add to.
 ** @param src The row to add.
 ** @param first Number of long integers to skip.
 ** @param len Number of long integers to add.
 ** @return Always returns dest.
 **/

PTR FfAddRowPartial(PTR dest, PTR src, int first, int len)
{
    register long i;

    if (FfChar == 2)	/* characteristic 2 is simple... */
#ifdef ASM_MMX
	__asm__("\n	movl 8(%ebp),%ecx\n"
		"	movl 12(%ebp),%ebx\n"
		"	movl 16(%ebp),%edx\n"
		"       sall $2,%edx\n"
		"       addl %edx,%ecx\n"
		"       addl %edx,%ebx\n"
		"       movl 20(%ebp),%edx\n"
		"	sarl $1,%edx\n"
		"	je .ADDROWPART_1\n"
		"	.align 16\n"
		".ADDROWPART_2:\n"
		"	movq (%ebx),%mm0\n"
		"	addl $8,%ebx\n"
		"	pxor (%ecx),%mm0\n"
		"	movq %mm0,(%ecx)\n"
		"	addl $8,%ecx\n"
		"	decl %edx\n"
		"	jne .ADDROWPART_2\n"
		".ADDROWPART_1:\n"
	       );
#else
    {	register long *l1 = (long *) dest + first;
	register long *l2 = (long *) src + first;
	for (i = len; i != 0; --i)
	{
	    register long x = *l2++;
	    *l1 ^= x;
	    l1++;
	}
    }
#endif
    else		/* any other characteristic */
    {	register BYTE *p1 = dest + first * sizeof(long);
	register BYTE *p2 = src + first * sizeof(long);
	for (i = len*sizeof(long); i != 0; --i)
	{
	    register int x = *p2++;
	    *p1 = mtx_tadd[*p1][x];
	    p1++;
	}
    }
    return dest;
}




/**
 ** Multiply a row by a coefficient.
 ** This function multiplies each element of @em row by @em mark.
 ** The row size and field order must have been set before.
 ** Multiplying a row with zero (FF_ZERO) initializes all elements to zero
 ** and is permitted even if @em row points into uninitialized memory.
 **/

void FfMulRow(PTR row, FEL mark)

{
    register BYTE *m;
    register BYTE *multab;
    register int i;

    CHECKFEL(mark);
    if (mark == FF_ZERO)
    {
        register long *l = (long *)row;
        for (i = LPR; i > 0; --i) *l++ = 0;
    }
    else if (mark != FF_ONE)
    {
	multab = mtx_tmult[mark];
	m = row;
	for (i = FfTrueRowSize(FfNoc); i != 0; --i)
	{
	    register int x = *m;
	    if (x != 0) *m = multab[x];
	    ++m;
	}
    }
}



/**
 ** Add a multiple of a row.
 ** This function adds a multiple of @em src to @em dest.
 **/

void FfAddMulRow(PTR dest, PTR src, FEL f)
{
    register int i;
    register BYTE *p1, *p2, *multab;

    CHECKFEL(f);
    if (f == FF_ZERO)
	return;
    if (f == FF_ONE)
    {
	FfAddRow(dest,src);
	return;
    }
    multab = mtx_tmult[f];
    p1 = dest;
    p2 = src;
    for (i = FfTrueRowSize(FfNoc); i != 0; --i)
    {
	*p1 = mtx_tadd[*p1][multab[*p2++]];
	++p1;
    }
}


/**
 ** Convert integer to field element.
 ** This function converts an integer to a field element using the same mapping as explained
 ** under FfToInt(). Both functions are inverse in the sense that the expression
 ** <tt>f == FfFromInt(FfToInt(f))</tt> is always true for valid field elements.
 ** FfFromInt() should be used whenever field elements are to be read with scanf().
 **/

FEL FfFromInt(int l)
{
    register int f;
    f = l % FfOrder;
    if (f < 0) f += FfOrder;
    return (FEL) f;
}

/**
 ** Convert field element to integer.
 ** This function converts a field element to an integer, using a "canonical" representation
 ** of field elements as integers which may be different from the internal representation.
 ** In particular, the prime field is mapped on {0,...p-1} with 0 representing the zero element
 ** and 1 the unit element. FfToInt() should be used whenever field elements are to be written
 ** with printf().
 **/

int FfToInt(FEL f)
{	
    return (int) f;
}


/**
 ** Multiply a vector by a matrix.
 ** This function multiplies the vector @em row from the right by the matrix @em mat and
 ** stores the result into @em result.
 ** The number of columns in both @em mat and @em result is determined by the current row size.
 ** (see FfNoc()).
 ** @attention @em result and @em row must not overlap. Otherwise the result is
 ** undefined.
 ** @param row The source vector (FfNoc columns).
 ** @param matrix The matrix (nor by FfNoc).
 ** @param nor number of rows in the matrix.
 ** @param[out] result The resulting vector (@em nor columns).
 **/

void FfMapRow(PTR row, PTR matrix, int nor, PTR result)

{
    register int i;
    register FEL f;
    BYTE *m = (BYTE *) matrix;
    register long *l = (long *)result;

#ifdef DEBUG
    if (result >= row && result < row + FfRowSize(nor))
	MTX_ERROR("row and result overlap: undefined result!");
    if (row >= result && row < result + FfCurrentRowSize)
	MTX_ERROR("row and result overlap: undefined result!");
#endif

    /* Fill the result with zeroes.
       ---------------------------- */
    for (i = LPR; i > 0; --i)
	*l++ = 0;

    if (FfOrder == 2)       /* GF(2) is a special case */
    {
        register long *x1 = (long *) matrix;
        register BYTE *r = (BYTE *) row;

        for (i = nor; i > 0; ++r)
        {
	    register BYTE mask;
	    if (*r == 0)
	    {
		i -= 8;
		x1 += 8 * LPR;
		continue;
	    }
	    for (mask = 0x80; mask != 0 && i > 0; mask >>= 1, --i)
	    {
                if ((mask & *r) == 0)
		{
                    x1 += LPR;  /* Skip that row */
		    continue;
		}
#ifdef ASM_MMX
__asm__("    pushl %ebx\n");
__asm__("    movl %0,%%ebx" : : "g" (x1) );
__asm__("    pushl %ecx\n"
	"    pushl %edx\n"
	"    movl 20(%ebp),%ecx\n"	/* result */
	);
__asm__ (
        "    movl LPR,%edx\n"
        "    sarl $1,%edx\n"
        "    je .FASTXOR_1\n"
        "    .align 16\n"
	".FASTXOR_2:\n"
        "    movq (%ebx),%mm0\n"
        "    addl $8,%ebx\n"
        "    pxor (%ecx),%mm0\n"
        "    movq %mm0,(%ecx)\n"
        "    addl $8,%ecx\n"
        "    decl %edx\n"
        "    jne .FASTXOR_2\n"
	".FASTXOR_1:\n"
	"    popl %edx\n"
	"    popl %ecx\n");
__asm__("    movl %%ebx,%0" : : "g" (x1) );
__asm__("    popl %ebx\n"
	);

#else
		{
                register long *x2 = (long *)result;
                register int k;
                for (k = LPR; k; --k)
                    *x2++ ^= *x1++;
		}
#endif
            }
        }
    }
    else                /* Any other field */
    {
        register BYTE *brow = (BYTE *) row;
        register int pos = 0;

        for (i = nor; i > 0; --i)
        {
            f = mtx_textract[pos][*brow];
            if (++pos == (int) MPB)
            {
                pos = 0;
                ++brow;
            }
            if (f != FF_ZERO)
            {
                register BYTE *v = m;
                register BYTE *r = result;
                register int k = FfCurrentRowSizeIo;
                if (f == FF_ONE)
                {
                    for (; k != 0; --k)
                    {
                        *r = mtx_tadd[*r][*v++];
                        ++r;
                    }
                }
                else
                {
                    register BYTE *multab = mtx_tmult[f];
                    for (; k != 0; --k)
                    {
		 	if (*v != 0)
			    *r = mtx_tadd[multab[*v]][*r];
			++v;
                        ++r;
                    }
                }
            }
            m += FfCurrentRowSize;              /* next row */
        }
    }
}






/**
 ** Scalar Product of Two Vectors.
 ** Given two vectors @f$a=(a_i)@f$ and @f$b=(b_i)@f$, this function calculates the
 ** scalar product @f$p=\sum_ia_ib_i@f$.
 ** @param a The first vector.
 ** @param b The second vector.
 ** @return Scalar product of the two vectors.
 **/


FEL FfScalarProduct(PTR a, PTR b)
{
    register unsigned char *ap = (unsigned char *) a;
    register unsigned char *bp = (unsigned char *) b;
    register int i;
    FEL f = FF_ZERO;

    for (i = FfNoc-1; i >= MPB; i -= MPB)
    {

	register int k;
	if (*ap != 0 && *bp != 0)
	{
	    for (k = 0; k < MPB; ++k)
		f = FfAdd(f,FfMul(mtx_textract[k][*ap],mtx_textract[k][*bp]));
	}
	++ap;
	++bp;
    }
    for (; i >= 0; --i)
	f = FfAdd(f,FfMul(mtx_textract[i][*ap],mtx_textract[i][*bp]));

    return f;
}




/**
!section kernel.ff.row
 ** Extract one column of a matrix.
 ** @param mat
    Pointer to the matrix.
 ** @param nor
    Number of rows in matrix.
 ** @param col
    Column to extract (starting with 1).
 ** @param result
    Pointer to buffer for the extracted column.
!description
    This function extracts one column out of a matrix and stores it as
    a row vector in |result|. The number of columns of the matrix must have
    been set with |FfSetNoc()|. |nor| is the number of rows in the matrix. The
    result is a row with |nor| entries, i.e., the length of |result| must be at
    least |nor|. |mat| and |result| must not overlap, or the result is
    undefined.
 **/

void FfExtractColumn(PTR mat,int nor,int col,PTR result)

{
    register BYTE *x = (BYTE *)mat + (col / MPB);
    register BYTE *extab = mtx_textract[col % MPB];
    register BYTE a = 0;
    register int ind = 0;
    register BYTE *y = result;
    register int count;

    for (count = nor; count > 0; --count)
    {
	a = (BYTE) (a + mtx_tinsert[ind][extab[*x]]);
	if (++ind == MPB)
	{
	    *y++ = a;
	    a = 0;
	    ind = 0;
	}
	x += FfCurrentRowSize;
    }
    if (ind != 0) *y = a;
}





/**
 ** Insert a mark into a row
 ** This function inserts the field element @em mark at position @em col into @em row.
 ** Column indexes start with 0.
 ** Before this function can be used, the field must be selected with FfSetField().
 ** %FfInsert() does not need the row size beeing set correctly. For
 ** example, if you are working with rows of different size, you do
 ** not have to call FfSetNoc() prior to each %FfInsert(). On the
 ** other hand, there is no protection against writing beyond the end
 ** of a row.
 **
 ** If the MeatAxe is compiled with the DEBUG option
 ** %FfInsert() checks that @em mark is a valid
 ** field element and @em col is not negative. If also the PARANOID option
 ** was in effect during compilation, %FfInsert() also checks if @em col is
 ** less than or equal to the current row size.
 ** @param row Pointer to the row.
 ** @param col Insert position (0-based).
 ** @param mark Value to insert.
 **/

void FfInsert(PTR row, int col, FEL mark)
{
    register BYTE *loc = (BYTE *)row + (col / MPB);
    register int idx = col % MPB;

#ifdef DEBUG
    CHECKFEL(mark);
    if (col < 0)
	MTX_ERROR1("Invalid column index %d",col);
#ifdef PARANOID
    if (col >= FfNoc)
	MTX_ERROR2("Invalid column index %d, noc=%d)",col,FfNoc);
#endif
#endif
    *loc = (BYTE) (mtx_tnull[idx][*loc] + mtx_tinsert[idx][mark]);
}



/**
!function FfExtract  "Extract a mark from a row"
 ** @param row
    Pointer to the row.
 ** @param col
    Index of mark to extract (0-based).
 ** @return
    |col|-th entry of |row|.
!description
    This function returns the entry at position |col| of a
    row.  Note that column indexes start with 0, i.e., |FfExtract(row,0)|
    returns the first entry of a row. Like |FfInsert()|, this function
    does not depend on the current row size. Reading beyond the end of
    a row will probably not produce an error, but the result is undefined.
 ** @see FfInsert
 **/

FEL FfExtract(PTR row, int col)

{
#ifdef DEBUG
    if (col < 0)
	MTX_ERROR1("Invalid column index (%d)",col);
#ifdef PARANOID
    if (col >= FfNoc)
	MTX_ERROR2("Invalid column index %d, noc=%d",col,FfNoc);
    CHECKFEL(mtx_textract[col % MPB][((BYTE *)row)[col / MPB]]);
#endif
#endif
    return mtx_textract[col % MPB][((BYTE *)row)[col / MPB]];
}



/**
 ** Find pivot column.
 ** This function scans the vector @a row and finds the first non-zero
 ** mark. The mark is stored into <tt>*mark</tt> and its position (counting
 ** from 0) is returned. If the whole vector is zero, %FfFindPivot()
 ** returns -1 and leaves <tt>*mark</tt> unchanged.
 ** @param row Pointer to the row.
 ** @param mark Buffer for pivot element.
 ** @return Index of the first non-zero entry in @a row or -1 if all entries are zero.
 **/

int FfFindPivot(PTR row, FEL *mark)
{
    register long *l = (long *) row;
    register int idx;
    register BYTE *m;

    for (idx = 0; idx < LPR && *l == 0; ++idx, ++l);
    if (idx == LPR)
	return -1;
    idx = idx * sizeof(long) * MPB;
    m = (BYTE *)l;
    while (*m == 0)
    {	++m;
	idx += MPB;
    }
    idx += mtx_tffirst[*m][1];
    if (idx >= FfNoc)		/* Ignore garbage in padding bytes */
	return -1;
    *mark = mtx_tffirst[*m][0];
    return idx;
}


/**
 ** @}
 **/


