////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - library interfaces
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MEATAXE_H_INCLUDED
#define MEATAXE_H_INCLUDED

#include <stdio.h>
#include <stdarg.h>

// Version nameng convention:
// x.y.z          - Released version
// x.y.z-UNSTABLE - Future version under development
#define MTX_VERSION "2.5.0-UNSTABLE"


#if defined GCC
#define MTX_PRINTF_ATTRIBUTE(f,v) __attribute__((format(printf,f,v)))
#else
#define MTX_PRINTF_ATTRIBUTE(f,v)
#endif

extern char MtxVersion[];       /**< The MeatAxe version. */

/** @addtogroup os
 * @{
 **/

int sysCreateDirectory(const char *name);
int sysGetPid();
void sysInit(void);
void *sysMalloc(size_t nbytes);
FILE *sysFopen(const char *name, const char*mode);
void sysFree(void *x);
int sysFseek(FILE *f,long pos);
int sysFseekRelative(FILE *file, long distance);
void *sysRealloc(void *buf, size_t nbytes);
int sysReadLong32(FILE *f, long *buf, int n);
int sysReadLongX(FILE *f, long *buf, int n);
int sysRemoveDirectory(const char *name);
int sysRemoveFile(const char *name);
void sysSetTimeLimit(long nsecs);
long sysTimeUsed(void);
int sysWriteLong32(FILE *f, const long *buf, int n);
int sysWriteLongX(FILE *f, const long *buf, int n);

#define ALLOC(type) ((type *) sysMalloc(sizeof(type)))
#define NALLOC(type,n) ((type *) sysMalloc((size_t)(n) * sizeof(type)))
#define NREALLOC(x,type,n) \
   ((type *) sysRealloc(x,(size_t)(n) * sizeof(type)))
#define FREE(x) sysFree(x)

/** @} **/

/* ---------------------------------------------------------------------------------------------- */

/**
** @addtogroup ff
** @{
**/

/* Data types and constants
   ------------------------ */

#if MTXZZZ == 0

typedef unsigned char FEL;              /**< A finite field element */
typedef FEL *PTR;                       /**< A pointer to a row vector */
#define FF_ZERO ((FEL)0)                /**< The zero field element */
#define FF_ONE ((FEL)1)                 /**< The unit element */
#define ZZZVERSION 6

#elif MTXZZZ == 1

typedef unsigned short FEL;
typedef unsigned short *PTR;
#define FF_ZERO ((FEL)0xFFFF)
#define FF_ONE ((FEL)0)
#define ZZZVERSION 0x105

#else

#error "MTXZZZ undefined"

#endif

/* Kernel variables and functions
   ------------------------------ */

extern int ffOrder;
extern int ffChar;
extern FEL ffGen;
extern int ffNoc; // TODO: REMOVE

/* Arithmetic */
FEL ffAdd(FEL a, FEL b);
FEL ffSub(FEL a, FEL b);
FEL ffMul(FEL a, FEL b);
FEL ffDiv(FEL a, FEL b);
FEL ffNeg(FEL a);
FEL ffInv(FEL a);

int ffMakeTables(int field);
int ffSetField(int field);

//TODO: REMOVE
void ffSetNoc(int noc);

/// Add a multiple of a row.
/// This function adds a multiple of @em src to @em dest.
void ffAddMulRow(PTR dest, PTR src, FEL f);

void ffAddMulRowPartial(PTR dest, PTR src, FEL f, int firstcol);

/// Add two rows.
/// This function adds src to dest. Field order and row size must have been set before.
/// @param dest The row to add to.
/// @param src The row to add.
/// @return Always returns dest.
PTR ffAddRow(PTR dest, PTR src);
PTR ffAddRowPartial(PTR dest, PTR src, int first);

PTR ffAlloc(int nor, int noc);

int ffCmpRows(PTR p1, PTR p2);
void ffCleanRow(PTR row, PTR matrix, int nor, int noc, const int *piv);
int ffCleanRow2(PTR row, PTR matrix, int nor, int noc, const int *piv, PTR row2);

/// Clean Row and Repeat Operations.
/// This function works like ffCleanRow(), but repeats all operations on
/// a second row/matrix.
/// @param row Pointer to row to be cleaned.
/// @param mat Matrix to clean with.
/// @param nor Number of rows.
/// @param piv Pivot table for @em mat.
/// @param row2 Pointer to the second row to be cleaned.
/// @param mat2 Matrix to the second matrix.
/// @return 0 on success, -1 on error.
int ffCleanRowAndRepeat(PTR row, PTR mat, int nor, int noc, const int *piv, PTR row2, PTR mat2);

void ffCopyRow(PTR dest, PTR src);

/// Embed a subfield.
/// @param a Element of the subfield.
/// @param subfield Subfield order. Must be a divisor of the current field order.
/// @return @em a, embedded into the current field. (FEL)255 on error.
FEL ffEmbed(FEL a, int subfield);

FEL ffExtract(PTR row, int col);
void ffExtractColumn(PTR mat,int nor,int col,PTR result);
int ffFindPivot(PTR row, FEL *mark);
void ffFree(PTR x);

/// Convert integer to field element.
/// This function, together with ffFromInt(), defines a bijection between field elements and
/// the set of integers {0, 1, ... q-1}, where q is the field order. This mapping is used when
/// field elements are converted to or from text form and has the following properties:
/// - ffFromInt(0) if the zero element
/// - ffFromInt(1) if the unit element
/// - The restriction to {0,..,p-1} (with the usual arithmetic mod p) is an isomorphism of
///   Z/(pZ) to the prime field.
FEL ffFromInt(int l);

/// Get row pointer.
/// This function returns a pointer to a row of a matrix, given the row index.
/// @a base must be a pointer to the beginning of a row, but this need not be the first
/// row of the matrix. For example, <tt>x = ffGetPtr(x,1,noc)</tt> can be used to advance a
/// row pointer to the next row.
///
/// Note: The function does not check if the resulting pointer is still inside the matrix.
/// @see ffStepPtr()
/// @param base Pointer to the first row of the matrix.
/// @param row Row index. The first row has index 0.
PTR ffGetPtr(PTR base, int row, int noc);

void ffMulRow(PTR row, FEL mark);

FILE *ffReadHeader(const char *name, int *fld, int *nor, int *noc);
int ffReadRows(FILE *f, PTR buf, int n, int noc);
FEL ffRestrict(FEL a, int subfield);
ssize_t ffSize(int nor, int noc);
size_t ffRowSize(int noc);
FEL ffScalarProduct(PTR a, PTR b);
int ffSeekRow(FILE *f, int pos);
void ffStepPtr(PTR *x, int noc);
void ffSwapRows(PTR dest, PTR src);
const char *ffToGap(FEL f);

/// Convert field element to integer. See @ref FfFromInt for more information.
int ffToInt(FEL f);

size_t ffTrueRowSize(int noc);
FILE *ffWriteHeader(const char *name, int fld, int nor, int noc);
int ffWriteRows(FILE *f, PTR buf, int n, int noc);

/* --------------------------------------------------------------------------
   Macro versions of kernel functions
   -------------------------------------------------------------------------- */

#if MTXZZZ == 0

extern FEL mtx_tmult[256][256];
extern FEL mtx_tadd[256][256];
extern FEL mtx_taddinv[256], mtx_tmultinv[256];
extern FEL mtx_tffirst[256][2];
extern FEL mtx_textract[8][256];
extern FEL mtx_tnull[8][256];
extern FEL mtx_tinsert[8][256];
extern long mtx_embedord[4];
extern FEL mtx_embed[4][16];
extern FEL mtx_restrict[4][256];

#define ffAdd(a,b) ((FEL)mtx_tadd[(int)(unsigned char)a][(int)(unsigned char)b])
#define ffDiv(a,b) ffMul((a),ffInv(b))
#define ffInv(a) (mtx_tmultinv[(int)(unsigned char)a])
#define ffMul(a,b) ((FEL)mtx_tmult[(int)(unsigned char)a][(int)(unsigned char)b])
#define ffNeg(a) (mtx_taddinv[(int)(unsigned char)a])
#define ffSub(a,b) ffAdd((a),ffNeg(b))

void ffInsert(PTR row, int col, FEL mark);

#elif MTXZZZ == 1

#define ffExtract(row,col) ((FEL)((row)[col]))
#define ffInsert(row,col,mark) ((void)((row)[col] = mark))

#else
   #error Illegal value of MTXZZZ
#endif

/** @} **/

/* ------------------------------------------------------------------
   Other low-level functions (zzz2.c)
   ------------------------------------------------------------------ */

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Multiply a vector by a matrix.
/// This function multiplies the vector @a row from the right by the matrix @a mat and
/// stores the result into @a result.
/// The number of columns in both @a mat and @a result must be the current row size
/// (see @ref ffNoc). The number of rows in @a mat and the number of columns in @a row must be
/// equal to @a nor.
///
/// @attention @em result and @em row must not overlap. Otherwise the result is undefined.
///
/// @param row The source vector (ffNoc columns).
/// @param matrix The matrix (nor by ffNoc).
/// @param nor number of rows in the matrix.
/// @param[out] result The resulting vector (@em nor columns).
void ffMapRow(PTR row, PTR matrix, int nor, PTR result);

void ffPermRow(PTR row, const long *perm, PTR result);
int ffSumAndIntersection(int noc, PTR wrk1, int *nor1, int *nor2, PTR wrk2, int *piv);


////////////////////////////////////////////////////////////////////////////////////////////////////
// Library initialization and cleanup
////////////////////////////////////////////////////////////////////////////////////////////////////

extern int Mtx_IsInitialized;
extern int Mtx_IsBigEndian;
extern int MtxOpt_UseOldWordGenerator;

/// This function initializes the MeatAxe library including finite field arithmetic and file i/o
/// functions. It must be called before any other MeatAxe library function.
/// It is legal to call MtxInitLibrary() multiple times. Only the first call will actually do
/// anything.
/// An application that uses @ref MtxInitApplication need not call this function.
///
/// @a argv0 is the name of the process executable. It will be used to initialize directory names
/// such as @ref MtxLibDir, which have a default value relative to the executable directory. If the
/// program name is not known, the argument may be NULL or an empty string.
///  
/// @c MtxInitLibrary() returns a version number which is different for each implementation of the
/// arithmetic, or -1 on error.
int MtxInitLibrary(char* argv0);
void MtxCleanupLibrary();

/// Sets the library directory.
void mtxSetLibraryDirectory(const char *dir);

/// This variable contains the name of the MeatAxe library directory.
/// MtxLibDir can be set on the command line with the "-L" option. Otherwise, the value of the 
/// environment variable MTXLIB is used. If neither "-L" nor MTXLIB are defined, the directory is
/// assumed to be on the same level as the program execuable directory and named "lib". For example,
/// if the program is "/home/user1/mtx/bin/zcp", the derived library directory would be
/// "/home/user1/mtx/lib".
/// If this fails, too, the current directory (".") is used.
extern char MtxLibDir[];

/// @addtogroup str
/// @{

/// A dynamic string.
/// The only member, S, is a normal «char*» pointing to a NUL terminated text.
/// Note however, that dynamic strings use their own memory management
/// which cannot be mixed with the standard C library memory functions.
/// Unused strings must be freed with strFree(), and you must never use
/// free() or realloc() on a dynamic string.

typedef struct {
   char *S;     /* pointer to NUL terminated string */
} String;

String strAlloc(size_t initial_capacity);
void strFree(String *s);
void strAppend(String *s, const char *text);
MTX_PRINTF_ATTRIBUTE(2,3)
void strAppendF(String *s, const char *fmt, ...);
MTX_PRINTF_ATTRIBUTE(2,3)
void strPrintF(String *s, const char *fmt, ...);

/// @}

////////////////////////////////////////////////////////////////////////////////////////////////////

#define APP_MAX_ARGS 50

/** Application information structure.
   This data structure is used to store information about the application.
   It is used by the command line parser, e.g., to display the help
   text.
   \n Name is the program name,
   \n Description is a one-line description of the program.
   \n Help is a help text that is to be displayed when
   the user invokes the program with '--help'.
   Here is an example:
   <PRE>
    MtxApplicationInfo_t AppInfo = { "sample", "MeatAxe sample",
        "\nSYNTAX\n"
        "    sample [-a] [-l <level>] <input> <ouput>\n"
        "\nOPTIONS\n"
        "    -a, --all ........ Output all data\n"
        "    -l, --level ...... Set output level (default: 42)\n"
     };
   </PRE>
   \sa AppAlloc
 **/

typedef struct {
   const char *Name;            /**< Program name. */
   const char *Description;     /**< One-line description of the program. */
   const char *Help;            /**< Help text. */
} MtxApplicationInfo_t;

/// Application data.
/// This data structure stores all internal data needed by the Application support functions,
/// such as command line arguments, temporary directory names, and more.
///
/// See also @ref AppAlloc

typedef struct {
   MtxApplicationInfo_t const *AppInfo;         /**< Program name and description. */
   int OrigArgC;                                /**< Original argc from main(). */
   char **OrigArgV;                             /**< Original argv from main(). */
   int ArgC;                                    /**< Number of arguments. */
   char **ArgV;                                 /**< Arguments. */
   int OptEnd;                                  /**< Used internally. */
   unsigned long IsDone[APP_MAX_ARGS];          /**< Used internally. */
   const char *OptArg;                          /**< Used internally. */
   int OptInd;                                  /**< Used internally. */
   char TempDirName[200];                       /**< Directory fr temporary files. */
} MtxApplication_t;

MtxApplication_t *appAlloc(MtxApplicationInfo_t const *ai, int argc, char **argv);
int appFree(MtxApplication_t *a);
int appGetOption(MtxApplication_t *app, const char *spec);
int appGetCountedOption(MtxApplication_t *app, const char *spec);
const char *appGetTextOption(MtxApplication_t *app, const char *spec,
                             const char *dflt);
int appGetIntOption(MtxApplication_t *app, const char *spec, int dflt,
                    int min, int max);
int appGetArguments(MtxApplication_t *app, int min_argc, int max_argc);
const char *appCreateTempDir(MtxApplication_t *app);

#define MTX_COMMON_OPTIONS_SYNTAX \
   "[<Options>]"

#define MTX_COMMON_OPTIONS_DESCRIPTION \
   "    -Q ...................... Quiet, no messages\n" \
   "    -V ...................... Verbose, more messages\n" \
   "    -T <MaxTime> ............ Set CPU time limit [s]\n" \
   "    --version ............... Show version information\n"

/* ------------------------------------------------------------------
   Messages and error handling
   ------------------------------------------------------------------ */

// Error messages
extern const char MTX_ERR_NOMEM[];
extern const char MTX_ERR_GAME_OVER[];
extern const char MTX_ERR_DIV0[];
extern const char MTX_ERR_FILEFMT[];
extern const char MTX_ERR_BADARG[];
extern const char MTX_ERR_RANGE[];
extern const char MTX_ERR_NOTECH[];
extern const char MTX_ERR_NOTSQUARE[];
extern const char MTX_ERR_INCOMPAT[];
extern const char MTX_ERR_BADUSAGE[];
extern const char MTX_ERR_OPTION[];
extern const char MTX_ERR_NARGS[];
extern const char MTX_ERR_NOTMATRIX[];
extern const char MTX_ERR_NOTPERM[];

/// Defines a source code location. Used for error messages.
struct MtxSourceLocation {
   const char* file;    ///< The source file name.
   int line;            ///< The line number.
   const char* func;    ///< The function name.
};

/// The current source code location.
#define MTX_HERE &(struct MtxSourceLocation) \
   {.file = __FILE__, .line = __LINE__, .func = __func__}

/// Run-time error information.
struct MtxErrorInfo {
   const struct MtxSourceLocation* source;
   const char *message;
};

typedef void MtxErrorHandler_t(const struct MtxErrorInfo *);

MTX_PRINTF_ATTRIBUTE(2,3)
void mtxAbort(const struct MtxSourceLocation* sl, const char *text, ...);

MtxErrorHandler_t *MtxSetErrorHandler(MtxErrorHandler_t *h);

#define MTX_ASSERT(e, retval) do { \
      if (!(e)) { \
         mtxAbort(MTX_HERE,"Assertion failed: %s",# e); \
         return retval; \
      } \
   } while (0)

#ifdef MTX_DEBUG
#define MTX_ASSERT_DEBUG(e,retval) MTX_ASSERT(e,retval)
#define MTX_FALSE_DEBUG(e) MTX_FALSE(e)
#else
#define MTX_ASSERT_DEBUG(e,retval)
#define MTX_FALSE_DEBUG(e)
#endif

/* ------------------------------------------------------------------
   Messages
   ------------------------------------------------------------------ */

extern int MtxMessageLevel;
#define MSG0 (MtxMessageLevel >= 0)
#define MSG1 (MtxMessageLevel >= 1)
#define MSG2 (MtxMessageLevel >= 2)
#define MSG3 (MtxMessageLevel >= 3)
#define MSG4 (MtxMessageLevel >= 4)
#define MESSAGE(level,args) \
   (MtxMessageLevel >= (level) ? (printf args, fflush(stdout), 1) : 0)

MTX_PRINTF_ATTRIBUTE(2,3)
void mtxMessage(int level, const char* msg, ...);

/* ------------------------------------------------------------------
   Miscellaneous
   ------------------------------------------------------------------ */

void mtxRandomInit(unsigned seed);
long int mtxRandom(void);
int mtxRandomInt(int max);
long gcd(long a, long b);
long lcm(long a, long b);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Structured text files (stfXXX.c)
////////////////////////////////////////////////////////////////////////////////////////////////////

/// Structured text file.
/// This structure is used for reading from and writing to structured text files.
 
typedef struct {
   FILE *File;          /**< The stream we're using */
   char *LineBuf;       /**< Buffers one 'line' */
   char *GetPtr;        /**< Current input position */
   int LineBufSize;     /**< Current buffer size */
   int OutPos;          /**< Number of chars in current line (writing only) */
   int LineNo;          /**< Current line number (reading and writing) */
} StfData;

int stfBeginEntry(StfData *f, const char *name);
int stfClose(StfData *f);
int stfEndEntry(StfData *f);
int stfGetInt(StfData *f, int *buf);
const char *stfGetName(StfData *f);
int stfGetString(StfData *f, char *buf, size_t bufsize);
int stfGetVector(StfData *f, int *bufsize, int *buf);
StfData *stfOpen(const char *name, const char* mode);
int stfPut(StfData *f, const char *text);
int stfPutInt(StfData *f, int value);
int stfPutString(StfData *f, const char *text);
int stfPutVector(StfData *f, int size, const int *value);
int stfReadLine(StfData *f);
int stfWriteValue(StfData *f, const char *name, const char *value);
int stfWriteInt(StfData *f, const char *name, int value);
int stfWriteString(StfData *f, const char *name, const char *value);
int stfWriteVector(StfData *f, const char *name, int size, const int *value);
int stfMatch(StfData *f, const char *pattern);

////////////////////////////////////////////////////////////////////////////////////////////////////
// MeatAxe binary data files
////////////////////////////////////////////////////////////////////////////////////////////////////

/// MeatAxe data file object.
typedef struct {
   unsigned long Magic;         /**< Used internally. */
   int Field;                   /**< Field order or type id. */
   int Nor;                     /**< Number of rows. */
   int Noc;                     /**< Number of columns. */
   FILE *File;                  /**< File to read frmo/write to. */
   char *Name;                  /**< File name. */
} MtxFile_t;

int mfIsValid(const MtxFile_t *file);
MtxFile_t *mfOpen(const char *name);
MtxFile_t *mfCreate(const char *name, int field, int nor, int noc);
int mfClose(MtxFile_t *file);
int mfReadLong(MtxFile_t *f, long *buf, int nrows);
int mfReadRows(MtxFile_t *f, PTR buf, int nrows);
int mfWriteLong(MtxFile_t *f, const long *buf, int count);
int mfWriteRows(MtxFile_t *f, PTR buf, int nrows);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Matrices over a finite field.
////////////////////////////////////////////////////////////////////////////////////////////////////

/// A matrix over a finite field.
typedef struct {
   unsigned long Magic;  ///< @private
   //void* next;           ///< @private
   //void** prev;          ///< @private
   int Field;            ///< Field order.
   int Nor;              ///< Number of rows.
   int Noc;              ///< Number of columns.
   PTR Data;             ///< Data, organized as array of rows.
   size_t RowSize;       ///< Size (in bytes) of one row.
   int *PivotTable;      ///< Pivot table (if matrix is in echelon form).
} Matrix_t;

Matrix_t *matAdd(Matrix_t *dest, const Matrix_t *src);
Matrix_t *matAddMul(Matrix_t *dest, const Matrix_t *src, FEL coeff);
Matrix_t *matAlloc(int field, int nor, int noc);
int matClean(Matrix_t *mat, const Matrix_t *sub);
int matCompare(const Matrix_t *a, const Matrix_t *b);
int matCopyRegion(Matrix_t *dest, int destrow, int destcol,
                  const Matrix_t *src, int row1, int col1, int nrows, int ncols);
Matrix_t *matCut(const Matrix_t *src, int row1, int col1, int nrows, int ncols);
Matrix_t *matCutRows(const Matrix_t *src, int row1, int nrows);
Matrix_t *matDup(const Matrix_t *src);
int matEchelonize(Matrix_t *mat);
int matFree(Matrix_t *mat);
PTR matGetPtr(const Matrix_t *mat, int row);
Matrix_t *matId(int fl, int nor);
Matrix_t *matInverse(const Matrix_t *src);
int matIsValid(const Matrix_t *m);
Matrix_t *matLoad(const char *fn);
Matrix_t *matMul(Matrix_t *dest, const Matrix_t *src);
Matrix_t *matMulScalar(Matrix_t *dest, FEL coeff);
int matNullity(const Matrix_t *mat);
int matNullity__(Matrix_t *mat);
Matrix_t *matNullSpace(const Matrix_t *mat);
Matrix_t *matNullSpace_(Matrix_t *mat, int flags);
Matrix_t *matNullSpace__(Matrix_t *mat);
int matOrder(const Matrix_t *mat);
int matPivotize(Matrix_t *mat);
Matrix_t *matPower(const Matrix_t *mat, long n);
void matPrint(const char *name, const Matrix_t *m);
Matrix_t *matRead(FILE *f);
int matSave(const Matrix_t *mat, const char *fn);
FEL matTrace(const Matrix_t *mat);
Matrix_t *matTransposed(const Matrix_t *src);
void matValidate(const struct MtxSourceLocation* sl, const Matrix_t *m);
int matWrite(const Matrix_t *mat, FILE *f);

/* For internal use only */
void mat_DeletePivotTable(Matrix_t *mat);

/* ------------------------------------------------------------------
   Greased matrices
   ------------------------------------------------------------------ */

/**
** Extraction table for greasing.
** This  structure is used internally for all greased matrix operations.
**/

typedef struct {
   long ***tabs;  /* tables for different remainders
                     of byte numbers mod grrows */
   int *nrvals;   /* number of values produced by each table */
   int nrtabs;    /* number of tables used */
} GrExtractionTable_t;

const GrExtractionTable_t *GrGetExtractionTable(int fl,int grrows);

/**
** A greased matrix.
**/

typedef struct {
   unsigned long Magic;
   int Field, Nor, Noc;
   int GrRows;                  /* Grease level (# of rows, 0=no grease) */
   int GrBlockSize;             /* Vectors per block (= Field^GrRows) */
   int NumVecs;                 /* Total number of vectors in <PrecalcData> */
   PTR PrecalcData;             /* Precalculated data */
   const GrExtractionTable_t
   *ExtrTab;                    /* Extraction table */
   int MPB;                     /* Number of marks per byte */
} GreasedMatrix_t;

void GrMapRow(PTR v,GreasedMatrix_t *M, PTR w);
GreasedMatrix_t *GrMatAlloc(const Matrix_t *m, int gr_rows);
int GrMatFree(GreasedMatrix_t *mat);
int GrMatIsValid(const GreasedMatrix_t *mat);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Permutations
////////////////////////////////////////////////////////////////////////////////////////////////////

/// A Permutation.
typedef struct {
   unsigned long Magic;   /**< Used internally. */
   int Degree;            /**< Degree of the permutation. */
   long *Data;            /**< Images of 0,1,2,... */
} Perm_t;

Perm_t *permAlloc(int deg);
int permCompare(const Perm_t *a, const Perm_t *b);
Perm_t *permDup(const Perm_t *src);
int permFree(Perm_t *p);
Perm_t *permInverse(const Perm_t *src);
int permIsValid(const Perm_t *p);
Perm_t *permLoad(const char *fn);
Perm_t *permMul(Perm_t *dest, const Perm_t *src);
int permOrder(const Perm_t *perm);
void permPrint(const char *name, const Perm_t *perm);
Perm_t *permPower(const Perm_t *p, int n);
Perm_t *permRead(FILE *f);
int permSave(const Perm_t *perm, const char *fn);
void permValidate(const struct MtxSourceLocation* src, const Perm_t *p);
int permWrite(const Perm_t *perm, FILE *f);

void Perm_ConvertOld(long *data, int len);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Polynomials over a finite field
////////////////////////////////////////////////////////////////////////////////////////////////////

/// A polynomial.
typedef struct {
   unsigned long Magic; /**< Used internally. */
   int Field;           /**< Field order. */
   int Degree;          /**< Degree of the polynomial. */
   FEL *Data;           /**< Coefficients. Degree+1 values, starting with the
                             constant term. */
   int BufSize;         /**< Used internally for memory management. */
}
Poly_t;

Poly_t *polAdd(Poly_t *dest, const Poly_t *src);
Poly_t *polAlloc(int fl, int n);
int polCompare(const Poly_t *a, const Poly_t *b);
Poly_t *polDerive(Poly_t *pol);
Poly_t *polDivMod(Poly_t *a, const Poly_t *b);
Poly_t *polDup(const Poly_t *p);
int polFree(Poly_t *p);
Poly_t *polGcd(const Poly_t *a, const Poly_t *b);
int polGcdEx(const Poly_t *a, const Poly_t *b, Poly_t **result);
int polIsValid(const Poly_t *p);
Poly_t *polMod(Poly_t *a, const Poly_t *b);
void Pol_Normalize(Poly_t *p);
Poly_t *polLoad(const char *fn);
Poly_t *polMul(Poly_t *dest, const Poly_t *src);
void polPrint(char *name, const Poly_t *p);
Poly_t *polRead(FILE *f);
int polSave(const Poly_t *pol, const char *fn);
void polValidate(const struct MtxSourceLocation* sl, const Poly_t *p);
int polWrite(const Poly_t *p, FILE *f);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Factored polynomials
////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct {
   unsigned long Magic; /**< Used internally. */
   int NFactors;        /**< Number of different irreducible factors. */
   int BufSize;         /**< Used internally for memory management. */
   Poly_t **Factor;     /**< List of irreducible factors. */
   int *Mult;           /**< Multiplicity of each factor. */
} FPoly_t;

FPoly_t *fpAlloc();
int fpFree(FPoly_t *x);
int fpIsValid(const FPoly_t *p);
FPoly_t *fpMul(FPoly_t *dest, const FPoly_t *src);
FPoly_t *fpMulP(FPoly_t *dest, const Poly_t *src, int pwr);
int fpPrint(const char *name, const FPoly_t *p);
void fpValidate(const struct MtxSourceLocation* src, const FPoly_t *p);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Bit strings
////////////////////////////////////////////////////////////////////////////////////////////////////

/// A bit string.
typedef struct {
   unsigned long Magic;   ///< Used internally.
   int Size;              ///< Number of bits. 
   int BufSize;           ///< Used internally for memory management. 
   long Data[1];          ///< The bits. The least significant bit comes first.
} BitString_t;

BitString_t *bsAlloc(int size);
int bsAnd(BitString_t *dest, const BitString_t *src);
int bsClear(BitString_t *bs, int i);
int bsClearAll(BitString_t *bs);
int bsCompare(const BitString_t *a, const BitString_t *b);
BitString_t *bsCopy(BitString_t *dest, const BitString_t *src);
BitString_t *bsDup(const BitString_t *src);
int bsFree(BitString_t *bs);
int bsIntersectionCount(const BitString_t *a, const BitString_t *b);
int bsIsSub(const BitString_t *a, const BitString_t *b);
int bsIsValid(const BitString_t *bs);
int bsMinus(BitString_t *dest, const BitString_t *src);
int bsOr(BitString_t *dest, const BitString_t *src);
void bsPrint(const char *name, const BitString_t *bs);
BitString_t *bsRead(FILE *f);
int bsSet(BitString_t* bs, int i);
int bsTest(const BitString_t *bs, int i);
void bsValidate(const struct MtxSourceLocation* src, const BitString_t *bs);
int bsWrite(BitString_t *bs, FILE *f);

#ifndef MTX_DEBUG

#define BS_BPL (sizeof(long) * 8)
#define bsSet(bs,i) ((bs)->Data[(i) / BS_BPL] |= 1L << ((i) % BS_BPL))
#define bsClear(bs,i) ((bs)->Data[(i) / BS_BPL] &= ~(1L << ((i) % BS_BPL)))
#define bsTest(bs,i) (((bs)->Data[(i) / BS_BPL] & (1L << ((i) % BS_BPL))) != 0 ? 1 : 0)

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Sets (subsets of Z).
////////////////////////////////////////////////////////////////////////////////////////////////////

/// A set of integers.
typedef struct {
   unsigned long Magic;         /**< Used internally. */
   int Size;                    /**< Number of elements. */
   int BufSize;                 /**< Used internally for memory management. */
   long *Data;                  /**< The elements in ascending order. */
} Set_t;

Set_t *setAlloc();
int setContains(const Set_t *set, long elem);
Set_t *setDup(const Set_t *s);
int setFree(Set_t *x);
int setInsert(Set_t *set, long elem);
int setIsValid(const Set_t *s);
int setPrint(char *name, const Set_t *s);
void setValidate(const struct MtxSourceLocation* src, const Set_t *s);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Integer matrices.
////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct {
   unsigned long Magic;
   int Nor;     ///< Number of rows
   int Noc;     ///< Number of colums
   long *Data;  ///< Marks (row by row)
} IntMatrix_t;

IntMatrix_t *imatAlloc(int nor, int noc);
int imatFree(IntMatrix_t *mat);
int imatIsValid(const IntMatrix_t *m);
IntMatrix_t *imatLoad(const char *fn);
IntMatrix_t *imatRead(FILE *f);
int imatSave(const IntMatrix_t *mat, const char *file_name);
int imatWrite(const IntMatrix_t *mat, FILE *f);

/* --------------------------------------------------------------------------
   Polymorphic objects
   -------------------------------------------------------------------------- */

void *XDup(void *a);
int XIsCompatible(void *a, void *b);
void XFree(void *a);
void *XInverse(void *a);
void *XLoad(const char *fn);
void XMul(void *a, void *b);
long XOrder(void *a);
void *XPower(void *a, int n);
int XSave(void *a, const char *fn);

/* --------------------------------------------------------------------------
   Matrix sets
   -------------------------------------------------------------------------- */

/**
** An element of a matrix set.
**/

typedef struct {
   Matrix_t *Matrix;
   int PivRow;
   int PivCol;
   FEL PivMark;
} MatrixSetElement_t;

/**
** A set of matrices.
**/
typedef struct {
   unsigned long Magic;
   int Len;
   MatrixSetElement_t *List;
} MatrixSet_t;

MatrixSet_t *msAlloc();
int msClean(const MatrixSet_t *set, Matrix_t *mat);
int msCleanAndAppend(MatrixSet_t *set, Matrix_t *mat);
int msFree(MatrixSet_t *set);
int msIsValid(const MatrixSet_t *set);

/* --------------------------------------------------------------------------
   Matrix representations
   -------------------------------------------------------------------------- */

typedef struct {
   unsigned long Magic;
   int NGen;
   Matrix_t **Gen;
} MatRep_t;

#define MR_COPY_GENERATORS  0x0001

int mrAddGenerator(MatRep_t *rep, Matrix_t *gen, int flags);
MatRep_t *mrAlloc(int ngen, Matrix_t **gen, int flags);
int mrChangeBasis(MatRep_t *rep, const Matrix_t *trans);
int mrIsValid(const MatRep_t *rep);
int mrFree(MatRep_t *rep);
MatRep_t *mrLoad(const char *basename, int ngen);
int mrSave(const MatRep_t *rep, const char *basename);
MatRep_t *mrTransposed(const MatRep_t *rep);

/* ------------------------------------------------------------------
   The word generator
   ------------------------------------------------------------------ */

typedef struct {
   const MatRep_t *Rep;         /**< The representation. **/
   Matrix_t *Basis[8];          /**< Products of the generators **/
   int N2[8];                   /**< Coefficients **/
   int *Description;            /**< Symbolic description of a word **/
} WgData_t;

WgData_t *wgAlloc(const MatRep_t *rep);
int *wgDescribeWord(WgData_t *b, long n);
int wgFree(WgData_t *b);
Matrix_t *wgMakeWord(WgData_t *b, long n);
void wgMakeFingerPrint(WgData_t *b, int fp[6]);
const char *wgSymbolicName(WgData_t *b, long n);

/* ------------------------------------------------------------------
   Spin-up, Split, Quotients, etc.
   ------------------------------------------------------------------ */

#define SF_FIRST        0x0001  /* Try only the first seed vector */
#define SF_EACH         0x0002  /* Try each seed vector */
#define SF_MAKE         0x0004  /* Try all 1-dimensional subspaces */
#define SF_SUB          0x0010  /* Try until finding a proper subspace */
#define SF_CYCLIC       0x0020  /* Try until finding a cyclic vector */
#define SF_COMBINE      0x0040  /* Combine the spans */
#define SF_SEED_MASK    0x000F
#define SF_MODE_MASK    0x00F0
#define SF_STD          0x0100  /* Spin up 'canonically' */

Matrix_t *QProjection(const Matrix_t *subspace, const Matrix_t *vectors);
Matrix_t *QAction(const Matrix_t *sub, const Matrix_t *gen);
Matrix_t *SAction(const Matrix_t *sub, const Matrix_t *gen);

typedef struct {
   int MaxSubspaceDimension;
   int MaxTries;
   int Result;
} SpinUpInfo_t;

int SpinUpInfoInit(SpinUpInfo_t *info);
Matrix_t *SpinUp(const Matrix_t *seed, const MatRep_t *rep, int flags,
                 IntMatrix_t **script, SpinUpInfo_t *info);
Matrix_t *SpinUpWithScript(const Matrix_t *seed, const MatRep_t *rep,
                           const IntMatrix_t *script);
int Split(Matrix_t *subspace, const MatRep_t *rep,
          MatRep_t **sub, MatRep_t **quot);

int ConvertSpinUpScript(IntMatrix_t *script);

Matrix_t *SpinUpWithPermutations(const Matrix_t *seed,
                                 int ngen,
                                 const Perm_t **gen,
                                 int flags,
                                 IntMatrix_t **script,
                                 SpinUpInfo_t *info);

/* ------------------------------------------------------------------
   Seed vector generator
   ------------------------------------------------------------------ */

long MakeSeedVector(const Matrix_t *basis, long lastno, PTR vec);

/* ------------------------------------------------------------------
   Miscellaneous algorithms
   ------------------------------------------------------------------ */

Matrix_t *matInsert_(Matrix_t *mat, const Poly_t *pol);
Matrix_t *matInsert(const Matrix_t *mat, const Poly_t *pol);
int IsSubspace(const Matrix_t *sub, const Matrix_t *space, int ngen);

Matrix_t *matTensor(const Matrix_t *m1, const Matrix_t *m2);
int MatrixToVector(const Matrix_t *mat, Matrix_t *vecs, int n);
Matrix_t *VectorToMatrix(Matrix_t *vecs, int n, int noc);
Matrix_t *TensorMap(Matrix_t *vec, const Matrix_t *a, const Matrix_t *b);

int StablePower(const Matrix_t *mat, int *pwr, Matrix_t **ker);
int StablePower_(Matrix_t *mat, int *pwr, Matrix_t **ker);

/* ------------------------------------------------------------------
   Polynomial factorization (Berlekamp algorithm)
   ------------------------------------------------------------------ */

FPoly_t *Factorization(const Poly_t *pol);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Characteristic and minimal polynomials
////////////////////////////////////////////////////////////////////////////////////////////////////

extern long CharPolSeed; // TODO: remove

Poly_t *charPolFactor(const Matrix_t *mat);
FPoly_t *charPol(const Matrix_t *mat);
Poly_t *minPolFactor(const Matrix_t *mat);
FPoly_t *minPol(const Matrix_t *mat);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Submodule lattice functions
////////////////////////////////////////////////////////////////////////////////////////////////////

/**
** @addtogroup cfinfo
** @{
**/

#define MAXGEN 20       /* Max. number of generators */
#define LAT_MAXCF 200   /* Max. number of composition factors */
#define MAXCYCL 30000   /* Max. number of cyclic submodules */
#define MAXDOTL 90000   /* Max. number of dotted lines */
#define MAXNSUB 20000   /* Max. number of submodules */
#define LAT_MAXBASENAME 100

typedef struct {
   long dim, num, mult;
   long idword;                 /* Identifying word */
   Poly_t *idpol;
   long peakword;               /* Peak word */
   Poly_t *peakpol;
   long nmount;                 /* Number of mountains */
   long ndotl;                  /* Number of dotted lines */
   long spl;                    /* Degree of splitting field */
}
CfInfo;

typedef struct {
   char BaseName[LAT_MAXBASENAME];      /**< Base name */
   int Field;                           /**< Field order */
   int NGen;                            /**< Number of generators */
   int NCf;                             /**< Number of irred. constituents */
   CfInfo Cf[LAT_MAXCF];                /**< Data for irred. constituents */
   int NSocles;                         /**< Loewy length */
   int *Socle;                          /**< Mult. of constituents in socles */
   int NHeads;                          /**< Number of radical layers */
   int *Head;                           /**< Mult. of constituents in Heads */
} Lat_Info;

int latReadInfo(Lat_Info *li, const char *basename);
int latWriteInfo(const Lat_Info *li);
const char *latCfName(const Lat_Info *li, int cf);
int latAddHead(Lat_Info *li, int *mult);
int latAddSocle(Lat_Info *li, int *mult);

#define LAT_RG_INVERT           0x0001  /* Invert generators */
#define LAT_RG_TRANSPOSE        0x0002  /* Transpose generators */
#define LAT_RG_STD              0x0004  /* Use standard form */

MatRep_t *latReadCfGens(Lat_Info *info, int cf, int flags);

/**
** @}
**/

/* ------------------------------------------------------------------
   Tensor condensation package
   ------------------------------------------------------------------ */

/**
** @ingroup tp
** Tensor condensation state.
** The %TkData_t structure is used by the tensor condensation algorithm to store
** state information. By use of tkReadInfo() and tkWriteInfo() a program can read
** the information from a file, and write back the data to a file if it was modified.
**/

typedef struct {
   char NameM[LAT_MAXBASENAME];         /**< Name of right factor */
   char NameN[LAT_MAXBASENAME];         /**< Name of left factor */
   int Dim;                             /**< Dimension of condensed module */
   int NCf;                             /**< Number of relevant constituents */
   int CfIndex[2][LAT_MAXCF];           /**< Constituent number */
} TkData_t;

int tkReadInfo(TkData_t *tki, const char *name);
int tkWriteInfo(TkData_t *tki, const char *name);

/* !!!!!!!!!!!!!!! 2.3 STUFF below !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* !!!!!!!!!!!!!!! 2.3 STUFF below !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* !!!!!!!!!!!!!!! 2.3 STUFF below !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/* ------------------------------------------------------------------
   Return codes
   ------------------------------------------------------------------ */

#define EXIT_OK         0       /* Exit code: normal end */
#define EXIT_ERR        1       /*            error */

/* ------------------------------------------------------------------
   Function operating on representations and vector spaces
   ------------------------------------------------------------------ */

int IsIsomorphic(const MatRep_t *rep1, const CfInfo *info1,
                 const MatRep_t *rep2, Matrix_t  **trans, int use_pw);
int MakeEndomorphisms(const MatRep_t *rep, const Matrix_t *nsp,
                      Matrix_t *endo[]);
Matrix_t *HomogeneousPart(MatRep_t *m, MatRep_t *s, Matrix_t *npw,
                          const IntMatrix_t *op, int dimends);

/* ------------------------------------------------------------------
   Lattice drawing functions
   ------------------------------------------------------------------ */

/**
** Node in a graphical lattice representation.
**/
typedef struct {
   double PosX, PosY;           /* Position [0..1] */
   unsigned long UserData;      /* User-defined attributes */
   int Layer;                   /* Layer number */
   double Score;                /* Used in optimization */
   int ScoreCount;
} LdNode_t;

/**
** Graphical lattice representation.
**/
typedef struct {
   int NNodes;
   LdNode_t *Nodes;
   int *IsSub;          /* Incidence relation, <NNodes> * <NNodes> entries */
   int *LayerNo;        /* Layer numbers */
   int NLayers;
} LdLattice_t;

#define LD_ISSUB(l,i,k) ((l)->IsSub[(i) * (l)->NNodes + (k)])

LdLattice_t *ldAlloc(int num_nodes);
int ldFree(LdLattice_t *l);
int ldAddIncidence(LdLattice_t *lat, int sub, int sup);
int ldSetPositions(LdLattice_t *l);

/* OLD STUFF */
int ChangeBasisOLD(const Matrix_t *M, int ngen, const Matrix_t *gen[],
                   Matrix_t *newgen[]);

#endif  /* !defined(_MEATAXE_H_) */

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
