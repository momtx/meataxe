////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - library interfaces
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MEATAXE_H_INCLUDED
#define MEATAXE_H_INCLUDED

#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#if defined(MTX_DEFAULT_THREADS)
#include <pthread.h>
#endif

// Version naming convention:
// x.y.z          - Released version
// x.y.z-UNSTABLE - Future version under development
#define MTX_VERSION "2.5.0-UNSTABLE"


#if defined __GNUC__
#define MTX_PRINTF(f,v) __attribute__((format(printf,f,v)))
#else
#define MTX_PRINTF(f,v)
#endif

enum MtxObjectType {
   MTX_TYPE_PERMUTATION = 0xFFFFFFFF,
   MTX_TYPE_POLYNOMIAL = 0xFFFFFFFE,
   MTX_TYPE_BITSTRING_FIXED = 0xFFFFFFFD,
   MTX_TYPE_BITSTRING_DYNAMIC = 0xFFFFFFFC,
   MTX_TYPE_WORD_GENERATOR = 0xFFFFFFFA,
   MTX_TYPE_MATREP = 0xFFFFFFF9,
   MTX_TYPE_INTMATRIX = 0xFFFFFFF8,
   MTX_TYPE_BINFILE = 0xFFFFFFF7,
   MTX_TYPE_STFILE = 0xFFFFFFF6,
   MTX_TYPE_CPSTATE = 0xFFFFFFF5,
   MTX_TYPE_STRBUF = 0xFFFFFFF4,
   MTX_TYPE_FPOLY = 0xFFFFFFF3,
   MTX_TYPE_LATINFO = 0xFFFFFFF2,
   MTX_TYPE_MATRIX = 0xFFFFFF01,
   MTX_TYPE_BEGIN = 0xFFFFFF00
};

////////////////////////////////////////////////////////////////////////////////////////////////////

/// @addtogroup mm
/// @{

void* mmAlloc(uint32_t typeId, size_t size);
void mmFree(void* obj, uint32_t typeId);
void mmLeakCheck();
uint32_t mmCheckpoint();
void mmRollback(uint32_t checkpoint);

/// @}

/// @addtogroup os
/// @{

int sysCreateDirectory(const char* name);
FILE* sysFopen(const char* name, const char*mode);
void sysFree(void* x);
int sysFseek(FILE *f, long pos);
int sysFseekRelative(FILE *file, long distance);
const char* sysGetExecutableName(const char* argv0);
int sysGetPid();
void sysInit(void);
void* sysMalloc(size_t nbytes);
size_t sysPad(size_t x, size_t unit);
void sysRead16(FILE *f, void* buf, size_t n);
void sysRead32(FILE *f, void* buf, size_t n);
void sysRead8(FILE *f, void* buf, size_t n);
void* sysRealloc(void* buf, size_t nbytes);
int sysRemoveDirectory(const char* name);
int sysRemoveFile(const char* name);
void sysSetTimeLimit(long nsecs);
uint64_t sysTime();
int sysTimeout(uint64_t* buf, unsigned intervalInSeconds);
long sysTimeUsed(void);
int sysTryRead32(FILE *f, void* buf, size_t n);
void sysWrite16(FILE *f, const void* buf, size_t n);
void sysWrite32(FILE *f, const void* buf, size_t n);
void sysWrite8(FILE *f, const void* buf, size_t n);

#define ALLOC(type) ((type *) sysMalloc(sizeof(type)))
#define NALLOC(type,n) ((type *) sysMalloc((size_t)(n) * sizeof(type)))
#define NREALLOC(x,type,n) \
   ((type *) sysRealloc(x,(size_t)(n) * sizeof(type)))

/// @}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Parallel execution (multithreading) support
////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct PexGroup PexGroup_t;

struct ErrorContextStack* pexContextStack();
PexGroup_t* pexCreateGroup();
void pexExecute(PexGroup_t* group, void (*f)(void* userData), void* userData);
void pexExecuteRange(PexGroup_t* group, void (*f)(void* userData, size_t begin, size_t end),
	void* userData, size_t begin, size_t end);
//void pexFinally(PexGroup_t* group, void(*f)(void* userData), void* userData);
void pexInit(int nThreads);
const char* pexLogPrefix();
int pexPoolSize();
MTX_PRINTF(1,2)
void pexSetThreadName(const char* name, ...);
int pexThreadNumber();
void pexThrottle(PexGroup_t* group, int* isEnabled, int loadFactor);
void pexShutdown();
void pexSleep(unsigned timeInMs);
unsigned pexThreadId();
void pexWaitAll();
void pexWait(PexGroup_t* group);


////////////////////////////////////////////////////////////////////////////////////////////////////
// Finite fields kernel - basic types
////////////////////////////////////////////////////////////////////////////////////////////////////

/// @addtogroup ff
/// @{

#if MTX_ZZZ == 0

typedef uint8_t FEL;            ///< A finite field element
typedef FEL *PTR;               ///< A pointer to a row vector
#define FF_ZERO ((FEL)0)        ///< The zero field element
#define FF_ONE ((FEL)1)         ///< The unit element
#define MTX_ZZZVERSION 6

#define MTX_MAXSUBFIELDORD 16	// Maximal order of subfields
#define MTX_MAXSUBFIELDS 4      // Maximal number of subfields

#elif MTX_ZZZ == 1

typedef uint16_t FEL;
typedef uint16_t* PTR;
#define FF_ZERO ((FEL)0xFFFF)
#define FF_ONE ((FEL)0)
#define MTX_ZZZVERSION 0x105

#define MTX_MAXSUBFIELDS 16     // Maximal number of subfields (14, actually)

#else

#error "MTX_ZZZ undefined"

#endif

/// @}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Messages and error handling
////////////////////////////////////////////////////////////////////////////////////////////////////

// Error messages
extern const char MTX_ERR_GAME_OVER[];
extern const char MTX_ERR_DIV0[];
extern const char MTX_ERR_FILEFMT[];
extern const char MTX_ERR_BADARG[];
extern const char MTX_ERR_RANGE[];
extern const char MTX_ERR_NOTECH[];
extern const char MTX_ERR_NOTSQUARE[];
extern const char MTX_ERR_INCOMPAT[];
extern const char MTX_ERR_OPTION[];
extern const char MTX_ERR_NOTMATRIX[];
extern const char MTX_ERR_NOTPERM[];

/// Describes a source code location. Used for error messages.
struct MtxSourceLocation {
   const char* file;    ///< The source file name.
   int line;            ///< The line number.
   const char* func;    ///< The function name.
};

/// The current source code location
#define MTX_HERE (&(const struct MtxSourceLocation){__FILE__, __LINE__, __func__})

/// Run-time error information.
struct MtxErrorInfo {
   struct MtxSourceLocation source;
   const char* message;
};

typedef void MtxErrorHandler_t(const struct MtxErrorInfo *);
typedef const char* MtxErrorContextProvider(void* userData);

/// @private
struct ErrorContext {
   struct MtxSourceLocation source;
   char* title;
   MtxErrorContextProvider* contextProvider;
   void* userData;
};

/// @private
struct ErrorContextStack {
   struct ErrorContext* stack;
   int capacity;
   int size;
};

MTX_PRINTF(2,3)
void mtxAbort(const struct MtxSourceLocation* sl, const char* text, ...);
MTX_PRINTF(2,3)
int mtxBegin(const struct MtxSourceLocation* sl, const char* s, ...);
int mtxBeginP(MtxErrorContextProvider ec, void* userData);
void mtxEnd(int id);

MtxErrorHandler_t* MtxSetErrorHandler(MtxErrorHandler_t* h);

#define MTX_ASSERT(e) do { \
      if (!(e)) { \
         mtxAbort(MTX_HERE,"Assertion failed: %s",# e); \
      } \
   } while (0)

#ifdef MTX_DEBUG
#define MTX_ASSERT_DEBUG(e) MTX_ASSERT(e)
#define MTX_FALSE_DEBUG(e) MTX_FALSE(e)
#else
#define MTX_ASSERT_DEBUG(e)
#define MTX_FALSE_DEBUG(e)
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
// Binary data files
////////////////////////////////////////////////////////////////////////////////////////////////////

/// Binary data file object.
typedef struct {
   void* next;           ///< @private
   void** prev;          ///< @private
   uint32_t seq;         ///< @private
   uint32_t typeId;      ///< @private
   uint32_t header[3]; ///< Last read/written object header.
   FILE *file;         ///< File handle.
   char* name;         ///< File name.
} MtxFile_t;

void mfClose(MtxFile_t* file);
MtxFile_t* mfCreate(const char* name, uint32_t field, uint32_t nor, uint32_t noc);
MtxFile_t* mfOpen(const char* name, const char* mode);
void mfRead32(MtxFile_t* file, void* buf, size_t nrows);
void mfRead8(MtxFile_t* file, void* buf, size_t nrows);
void mfSkip(MtxFile_t* file, size_t nBytes);
int mfTryRead32(MtxFile_t* file, void* buf, size_t nrows);
int mfTryReadHeader(MtxFile_t* file);
uint32_t mfReadHeader(MtxFile_t* file);
uint32_t mfObjectType(const MtxFile_t* file);
void mfValidate(const struct MtxSourceLocation* src, const MtxFile_t* file);
void mfWrite32(MtxFile_t* file, const void* buf, size_t count);
void mfWrite8(MtxFile_t* file, const void* buf, size_t count);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Finite fields kernel
////////////////////////////////////////////////////////////////////////////////////////////////////

/// @addtogroup ff
/// @{

extern uint32_t ffOrder;        // Current field order
extern int ffChar;              // Current field characteristic
extern FEL ffGen;               // Generator for the current field.

/// An invalid value. Used in places where a row/colum index is expected to signal that no value
/// is avalable.
#define MTX_NVAL ((uint32_t)0xFFFFFFFF)

// Arithmetic
FEL ffAdd(FEL a, FEL b);
FEL ffSub(FEL a, FEL b);
FEL ffMul(FEL a, FEL b);
FEL ffDiv(FEL a, FEL b);
FEL ffNeg(FEL a);
FEL ffInv(FEL a);

void ffAddMulRowPartial(PTR dest, PTR src, FEL f, uint32_t firstcol, uint32_t noc);
void ffAddMulRow(PTR dest, PTR src, FEL f, uint32_t noc);
void ffAddRowPartial(PTR dest, PTR src, uint32_t first, uint32_t noc);
PTR ffAddRow(PTR dest, PTR src, uint32_t noc);
PTR ffAlloc(int nor, int noc);
int ffCleanRow2(PTR row, PTR matrix, uint32_t nor, uint32_t noc, const uint32_t* piv, PTR row2);
int ffCleanRowAndRepeat(
      PTR row, PTR mat, int nor, int noc, const uint32_t* piv, PTR row2, PTR mat2);
void ffCleanRow(PTR row, PTR matrix, int nor, int noc, const uint32_t* piv);
int ffCmpRows(PTR p1, PTR p2, int noc);
void ffCopyRow(PTR dest, PTR src, int noc);
FEL ffEmbed(FEL a, int subfield);
void ffExtractColumn(PTR mat,int nor,int noc,int col,PTR result);
FEL ffExtract(PTR row, int col);
uint32_t ffFindPivot(PTR row, FEL *mark, int noc);
void ffFree(PTR x);
FEL ffFromInt(int l);
PTR ffGetPtr(PTR base, int row, int noc);
int ffMakeTables(int field);
void ffMulRow(PTR row, FEL mark, int noc);
void ffReadRows(MtxFile_t* file, PTR buf, uint32_t nor, uint32_t noc);
FEL ffRestrict(FEL a, int subfield);
size_t ffRowSize(uint32_t noc);
size_t ffRowSizeUsed(int noc);
FEL ffScalarProduct(PTR a, PTR b, int noc);
void ffSetField(int field);
ssize_t ffSize(uint32_t nor, uint32_t noc);
void ffStepPtr(PTR *x, int noc);
void ffSwapRows(PTR dest, PTR src, int noc);
int ffToInt(FEL f);
void ffWriteRows(MtxFile_t* file, PTR buf, uint32_t nor, uint32_t noc);

/// List of subfield orders, terminated with 0.
extern int mtx_subfields[17];

////////////////////////////////////////////////////////////////////////////////////////////////////
// Macro versions of kernel functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#if MTX_ZZZ == 0

extern FEL mtx_tmult[256][256];
extern FEL mtx_tadd[256][256];
extern FEL mtx_taddinv[256], mtx_tmultinv[256];
extern FEL mtx_tffirst[256][2];
extern FEL mtx_textract[8][256];
extern FEL mtx_tnull[8][256];
extern FEL mtx_tinsert[8][256];
extern FEL mtx_embed[4][16];
extern FEL mtx_restrict[4][256];

#define ffAdd(a,b) ((FEL)mtx_tadd[(uint8_t)a][(uint8_t)b])
#define ffDiv(a,b) ffMul((a),ffInv(b))
#define ffInv(a) (mtx_tmultinv[(uint8_t)a])
#define ffMul(a,b) ((FEL)mtx_tmult[(uint8_t)a][(uint8_t)b])
#define ffNeg(a) (mtx_taddinv[(uint8_t)a])
#define ffSub(a,b) ffAdd((a),ffNeg(b))

void ffInsert(PTR row, int col, FEL mark);

#elif MTX_ZZZ == 1

#define ffExtract(row,col) ((FEL)((row)[col]))
#define ffInsert(row,col,mark) ((void)((row)[col] = mark))

#else
   #error Illegal value of MTX_ZZZ
#endif

/// @}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Other low-level functions (zzz2.c)
////////////////////////////////////////////////////////////////////////////////////////////////////

void ffMapRow(PTR result, PTR row, PTR matrix, int nor, int noc);
void ffPermRow(PTR result, PTR row, const uint32_t* perm, int noc);
int ffSumAndIntersection(int noc, PTR wrk1, uint32_t* nor1, uint32_t* nor2, PTR wrk2, uint32_t* piv);


////////////////////////////////////////////////////////////////////////////////////////////////////
// Library initialization and cleanup
////////////////////////////////////////////////////////////////////////////////////////////////////

int mtxIsBigEndian();
void mtxCleanupLibrary();
void mtxInitLibrary(char* argv0);
const char* mtxLibraryDirectory();
const char* mtxVersion();

////////////////////////////////////////////////////////////////////////////////////////////////////

/// @addtogroup str
/// @{

/// A dynamic buffer used to construct strings.

typedef struct {
   void* next;           ///< @private
   void** prev;          ///< @private
   uint32_t seq;         ///< @private
   uint32_t typeId;      ///< @private
   size_t size;         // number of characters (not counting the terminating NUL)
   size_t capacity;     // max number of charcters (not counting the terminating NUL)
   char* data;
} StrBuffer_t;

StrBuffer_t* sbAlloc(size_t initialCapacity);
const char* sbData(StrBuffer_t* sb);
void sbAppend(StrBuffer_t* sb, const char* fragment);
void sbClear(StrBuffer_t* sb);
char* sbCopy(StrBuffer_t* sb);
MTX_PRINTF(1,2) char* strEprintf(const char* s, ...);
void sbFree(StrBuffer_t* sb);
MTX_PRINTF(2,3) void sbPrintf(StrBuffer_t* sb, const char* fmt, ...);
char* sbToString(StrBuffer_t* sb);
char* sbToEphemeralString(StrBuffer_t* sb);
void sbVprintf(StrBuffer_t* sb, const char* fmt, va_list args);

int strCompareRange(const char* x, const char* xEnd, const char* y, const char* yEnd);
char *strDup(const char* c);
const char* strPrefix(const char* s, const char* prefix);
char* strRange(const char* begin, const char* end);
char* strMakeEphemeral(char* c);
char* strVEprintf(const char* s, va_list args);
char* strVMprintf(const char* s, va_list args);
MTX_PRINTF(1,2) char* strMprintf(const char* s, ...);

/// @}

////////////////////////////////////////////////////////////////////////////////////////////////////

#define APP_MAX_ARGS 50

#define EXIT_OK         0       ///< Normal exit
#define EXIT_ERR        1       ///< Program aborted after error


/// Application information structure.
/// This data structure is used to store information about the application.
/// It is used by the command line parser, e.g., to display the help text.
/// See also @ref appAlloc

typedef struct {
   const char* name;            ///< Program name.
   const char* description;     ///< One-line description of the program.
   const char* help;            ///< Help text (shown with "--help").
} MtxApplicationInfo_t;

/// Application data.
/// This data structure stores all internal data needed by the Application support functions,
/// such as command line arguments, temporary directory names, and more.
///
/// See also @ref appAlloc

typedef struct {
   MtxApplicationInfo_t const* AppInfo; ///< Program name and description.
   int context;                         ///< Used internally
   int origArgC;                        ///< Original argc from main().
   char * const* origArgV;              ///< Original argv from main().
   int argC;                            ///< Number of arguments.
   char * const* argV;                  ///< Arguments.
   int optEnd;                          ///< Used internally.
   unsigned long isDone[APP_MAX_ARGS];  ///< Used internally.
   const char* optArg;                  ///< Used internally.
   char optName[100];                   ///< Used internally.
} MtxApplication_t;

MtxApplication_t* appAlloc(MtxApplicationInfo_t const* ai, int argc, char * const* argv);
void appFree(MtxApplication_t* a);
int appGetOption(MtxApplication_t* app, const char* spec);
const char* appGetTextOption(MtxApplication_t* app, const char* spec,
                             const char* dflt);
int appGetIntOption(MtxApplication_t* app, const char* spec, int dflt,
                    int min, int max);
int appGetArguments(MtxApplication_t* app, int min_argc, int max_argc);


#define MTX_COMMON_OPTIONS_SYNTAX \
   "[<Options>]"
#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)

#if defined(MTX_DEFAULT_THREADS)
   #define MTX_THREAD_OPTION_DESCRIPTION \
      "    -j <n> .................. Parallel execution on <n> CPU cores (default: "\
        STRINGIFY(MTX_DEFAULT_THREADS) ")\n"
#else
   #define MTX_THREAD_OPTION_DESCRIPTION \
      "    -j <n> .................. Ignored (threading support is disabled)\n"
#endif

#define MTX_COMMON_OPTIONS_DESCRIPTION \
   "    -Q ...................... Quiet, no messages\n" \
   "    -V ...................... Verbose, more messages\n" \
   "    -T <MaxTime> ............ Set CPU time limit [s]\n" \
   "    --log=[FILE]:LEVEL:[FMT]\n" \
   "                              Log to FILE (default: stdout) up to LEVEL (error,warning,\n" \
   "                              info,debug,debug2), using FORMAT (full,default).\n" \
   "    --help .................. Show help on command line syntax\n" \
   "    --version ............... Show version information\n" \
   MTX_THREAD_OPTION_DESCRIPTION

////////////////////////////////////////////////////////////////////////////////////////////////////
// Messages
////////////////////////////////////////////////////////////////////////////////////////////////////

enum {
   MTX_LOG_ERROR = -2,
   MTX_LOG_WARNING = -1,
   MTX_LOG_INFO = 0,
   MTX_LOG_DEBUG = 1,
   MTX_LOG_DEBUG2 = 2
};

/// @private
void logBuffered(StrBuffer_t* buf);

int logGetDefaultThreshold();
int logEnabled(int level);
void logInit(const char* spec);
void logPrepareForAbort();
MTX_PRINTF(2,3)
void logPrintf(int level, const char* msg, ...);
StrBuffer_t* logStart(int level);
void logSetDefaultThreshold(int level);

#define MTX_LOG(level, ...) \
   do { if (logEnabled(level)) { logPrintf((level), __VA_ARGS__);} } while(0)
#define MTX_LOGE(...) MTX_LOG(MTX_LOG_ERROR, __VA_ARGS__)
#define MTX_LOGW(...) MTX_LOG(MTX_LOG_WARNING, __VA_ARGS__)
#define MTX_LOGI(...) MTX_LOG(MTX_LOG_INFO, __VA_ARGS__)
#define MTX_LOGD(...) MTX_LOG(MTX_LOG_DEBUG, __VA_ARGS__)
#define MTX_LOG2(...) MTX_LOG(MTX_LOG_DEBUG2, __VA_ARGS__)

#define MTX_XLOG(level, sb) \
   for (StrBuffer_t* sb = logStart(level); sb; sb = (logBuffered(sb), NULL))
#define MTX_XLOGE(sb) MTX_XLOG(MTX_LOG_ERROR, sb)
#define MTX_XLOGW(sb) MTX_XLOG(MTX_LOG_WARNING, sb)
#define MTX_XLOGI(sb) MTX_XLOG(MTX_LOG_INFO, sb)
#define MTX_XLOGD(sb) MTX_XLOG(MTX_LOG_DEBUG, sb)
#define MTX_XLOG2(sb) MTX_XLOG(MTX_LOG_DEBUG2, sb)

////////////////////////////////////////////////////////////////////////////////////////////////////
// Miscellaneous
////////////////////////////////////////////////////////////////////////////////////////////////////

void mtxRandomInit(unsigned seed);
long int mtxRandom(void);
int mtxRandomInt(int max);
long gcd(long a, long b);
long lcm(long a, long b);
uint32_t gcd32u(uint32_t a, uint32_t b);
uint32_t lcm32u(uint32_t a, uint32_t b);
void hashLittle2(const void* data, size_t len, uint32_t* pc, uint32_t* pb);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Structured text files (stfXXX.c)
////////////////////////////////////////////////////////////////////////////////////////////////////

/// Structured text file.
/// This structure is used for reading from and writing to structured text files.

typedef struct {
   void* next;           ///< @private
   void** prev;          ///< @private
   uint32_t seq;         ///< @private
   uint32_t typeId;      ///< @private
   char* fileName;      ///< The file name.
   FILE *file;          ///< The stream we're using
   char* lineBuf;       ///< Buffers one 'line'
   char* getPtr;        ///< Current input position
   int lineBufSize;     ///< Current buffer size
   int outPos;          ///< Number of chars in current line (writing only)
   int lineNo;          ///< Current line number (reading and writing)
   int context;         ///< @private
} StfData;

int stfBeginEntry(StfData* f, const char* name);
void stfClose(StfData* f);
int stfEndEntry(StfData* f);
int stfGetInt(StfData* f, int* buf);
int stfGetU32(StfData* f, uint32_t* buf);
int stfGetULong(StfData* f, unsigned long* buf);
const char* stfGetName(StfData* f);
int stfGetString(StfData* f, char* buf, size_t bufsize);
int stfGetStringPtr(StfData* f, char** buf);
int stfGetVector(StfData* f, int* bufsize, int* buf);
StfData* stfOpen(const char* name, const char* mode);
int stfMatch(StfData* f, const char* pattern);
int stfPut(StfData* f, const char* text);
int stfPutInt(StfData* f, int value);
int stfPutString(StfData* f, const char* text);
int stfPutU32(StfData* f, uint32_t value);
int stfPutVector(StfData* f, int size, const int* value);
int stfReadLine(StfData* f);
void stfValidate(const struct MtxSourceLocation* where, const StfData* p);
int stfWriteInt(StfData* f, const char* name, int value);
int stfWriteString(StfData* f, const char* name, const char* value);
int stfWriteValue(StfData* f, const char* name, const char* value);
int stfWriteVector(StfData* f, const char* name, int size, const int* value);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Matrices over a finite field.
////////////////////////////////////////////////////////////////////////////////////////////////////

/// A matrix over a finite field.
typedef struct {
   void* next;           ///< @private
   void** prev;          ///< @private
   uint32_t seq;         ///< @private
   uint32_t typeId;      ///< @private
   uint32_t field;       ///< Field order.
   uint32_t nor;         ///< Number of rows.
   uint32_t noc;         ///< Number of columns.
   PTR data;             ///< Data, organized as array of rows.
   uint32_t* pivotTable; ///< Pivot table (if matrix is in echelon form).
} Matrix_t;

Matrix_t* matAdd(Matrix_t* dest, const Matrix_t* src);
Matrix_t* matAddMul(Matrix_t* dest, const Matrix_t* src, FEL coeff);
Matrix_t* matAlloc(int field, uint32_t nor, uint32_t noc);
uint32_t matClean(Matrix_t* mat, const Matrix_t* sub);
int matCompare(const Matrix_t* a, const Matrix_t* b);
void matCopyRegion(
   Matrix_t* dest, uint32_t destrow, uint32_t destcol,
   const Matrix_t* src, uint32_t row1, uint32_t col1, uint32_t nrows, uint32_t ncols);
Matrix_t* matCreateFromBuffer(PTR rows, int field, uint32_t nor, uint32_t noc);
Matrix_t* matDup(const Matrix_t* src);
Matrix_t* matDupRegion(
   const Matrix_t* src, uint32_t row0, uint32_t col0, uint32_t nrows, uint32_t ncols);
Matrix_t* matDupRows(const Matrix_t* src, uint32_t row1, uint32_t nrows);
uint32_t matEchelonize(Matrix_t* mat);
void matFree(Matrix_t* mat);
PTR matGetPtr(const Matrix_t* mat, uint32_t row);
Matrix_t* matId(int fl, uint32_t nor);
Matrix_t* matInverse(const Matrix_t* src);
int matIsValid(const Matrix_t* m);
Matrix_t* matLoad(const char* fn);
Matrix_t* matMul(Matrix_t* dest, const Matrix_t* src);
Matrix_t* matMulScalar(Matrix_t* dest, FEL coeff);
uint32_t matNullity(const Matrix_t* mat);
uint32_t matNullity__(Matrix_t* mat);
Matrix_t* matNullSpace(const Matrix_t* mat);
Matrix_t* matNullSpace_(Matrix_t* mat, int flags);
Matrix_t* matNullSpace__(Matrix_t* mat);
int matOrder(const Matrix_t* mat);
void matPivotize(Matrix_t* mat);
Matrix_t* matPower(const Matrix_t* mat, long n);
void matPrint(const char* name, const Matrix_t* m);
Matrix_t* matReadData(MtxFile_t* f);
void matSave(const Matrix_t* mat, const char* fn);
FEL matTrace(const Matrix_t* mat);
Matrix_t* matTransposed(const Matrix_t* src);
void matValidate(const struct MtxSourceLocation* sl, const Matrix_t* m);
void matWrite(const Matrix_t* mat, MtxFile_t* file);

/* For internal use only */
void mat_DeletePivotTable(Matrix_t* mat);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Permutations
////////////////////////////////////////////////////////////////////////////////////////////////////

/// A Permutation.
typedef struct {
   void* next;          ///< @private
   void** prev;         ///< @private
   uint32_t seq;        ///< @private
   uint32_t typeId;     ///< @private
   uint32_t degree;     ///< Degree of the permutation.
   uint32_t* data;      ///< Images of 0,1,2,...
} Perm_t;

Perm_t* permAlloc(uint32_t deg);
int permCompare(const Perm_t* a, const Perm_t* b);
void permConvertLegacyFormat(uint32_t* data, uint32_t degree);
Perm_t* permDup(const Perm_t* src);
void permFree(Perm_t* p);
Perm_t* permInverse(const Perm_t* src);
int permIsValid(const Perm_t* p);
Perm_t* permLoad(const char* fn);
Perm_t* permMul(Perm_t* dest, const Perm_t* src);
uint32_t permOrder(const Perm_t* perm);
Perm_t* permPower(const Perm_t* p, int n);
void permPrint(const char* name, const Perm_t* perm);
Perm_t* permRead(MtxFile_t* f);
Perm_t* permReadData(MtxFile_t* f);
void permSave(const Perm_t* perm, const char* fileName);
void permValidate(const struct MtxSourceLocation* src, const Perm_t* p);
void permWrite(const Perm_t* perm, MtxFile_t* file);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Polynomials over a finite field
////////////////////////////////////////////////////////////////////////////////////////////////////

/// A polynomial over a finite field.
typedef struct {
   void* next;          ///< @private
   void** prev;         ///< @private
   uint32_t seq;        ///< @private
   uint32_t typeId;     ///< @private
   uint32_t field;      ///< Field order.
   int32_t degree;      ///< Degree of the polynomial.
   FEL* data;           ///< Coefficients. Degree+1 values, starting with the constant term.
   uint32_t bufSize;    ///< Used internally for memory management.
}
Poly_t;

Poly_t* polAdd(Poly_t* dest, const Poly_t* src);
Poly_t* polAlloc(uint32_t field, int32_t degree);
int polCompare(const Poly_t* a, const Poly_t* b);
Poly_t* polDerive(Poly_t* pol);
Poly_t* polDivMod(Poly_t* a, const Poly_t* b);
Poly_t* polDup(const Poly_t* p);
void polFormat(StrBuffer_t* sb, const Poly_t* p);
void polFree(Poly_t* p);
Poly_t* polGcd(const Poly_t* a, const Poly_t* b);
int polGcdEx(const Poly_t* a, const Poly_t* b, Poly_t** result);
int polIsValid(const Poly_t* p);
Poly_t* polLoad(const char* fn);
Poly_t* polMod(Poly_t* a, const Poly_t* b);
Poly_t* polMul(Poly_t* dest, const Poly_t* src);
void polNormalize(Poly_t* p);
void polPrint(char* name, const Poly_t* p);
Poly_t* polReadData(MtxFile_t* f);
Poly_t* polRead(MtxFile_t* f);
void polSave(const Poly_t* pol, const char* fn);
char* polToEphemeralString(const Poly_t* p);
void polValidate(const struct MtxSourceLocation* sl, const Poly_t* p);
void polWrite(const Poly_t* p, MtxFile_t* file);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Factored polynomials
////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct {
   void* next;          ///< @private
   void** prev;         ///< @private
   uint32_t seq;        ///< @private
   uint32_t typeId;     ///< @private
   uint32_t field;
   uint32_t nFactors;   ///< Number of different irreducible factors.
   uint32_t bufSize;    ///< @private
   Poly_t** factor;     ///< List of irreducible factors.
   int* mult;           ///< Multiplicity of each factor.
} FPoly_t;

FPoly_t* fpAlloc(uint32_t field);
int fpCompare(const FPoly_t* a, const FPoly_t* b);
void fpFormat(StrBuffer_t* sb, const FPoly_t* p);
void fpFree(FPoly_t* x);
int fpIsValid(const FPoly_t* p);
FPoly_t* fpMul(FPoly_t* dest, const FPoly_t* src);
FPoly_t* fpMulP(FPoly_t* dest, const Poly_t* src, int pwr);
void fpPrint(const char* name, const FPoly_t* p);
char* fpToEphemeralString(const FPoly_t* p);
void fpValidate(const struct MtxSourceLocation* src, const FPoly_t* p);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Bit strings
////////////////////////////////////////////////////////////////////////////////////////////////////

/// A bit string.
typedef struct {
   void* next;           ///< @private
   void** prev;          ///< @private
   uint32_t seq;         ///< @private
   uint32_t typeId;      ///< @private
   size_t size;           ///< Number of signicant bits. Only used for fixed-size bit strings!
   size_t capacity;       ///< Maximum size.
   uint8_t* data;         ///< The bits. Bit 0 is the LSB of data[0].
} BitString_t;

BitString_t* bsAllocEmpty(void);
BitString_t* bsAlloc(size_t size);

void bsTrim(const BitString_t* bs);
void bsResize(BitString_t* bs, size_t newSize);

int bsFirst(const BitString_t* bs, size_t* i);
int bsNext(const BitString_t* bs, size_t* i);

void bsAnd(BitString_t* dest, const BitString_t* src);
void bsClear(BitString_t* bs, size_t i);
void bsClearAll(BitString_t* bs);
int bsCompare(const BitString_t* a, const BitString_t* b);
void bsCopy(BitString_t* dest, const BitString_t* src);
BitString_t* bsDup(const BitString_t* src);
int bsFree(BitString_t* bs);
size_t bsIntersectionCount(const BitString_t* a, const BitString_t* b);
int bsIsSub(const BitString_t* a, const BitString_t* b);
int bsIsValid(const BitString_t* bs);
void bsMinus(BitString_t* dest, const BitString_t* src);
void bsOr(BitString_t* dest, const BitString_t* src);
void bsPrint(const char* name, const BitString_t* bs);
BitString_t* bsReadData(MtxFile_t* f);
BitString_t* bsRead(MtxFile_t* f);
void bsSet(BitString_t* bs, size_t i);
void bsSkip(MtxFile_t* f);
int bsTest(const BitString_t* bs, size_t i);
void bsValidate(const struct MtxSourceLocation* src, const BitString_t* bs);
void bsWrite(const BitString_t* bs, MtxFile_t* file);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Integer matrices.
////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct {
   void* next;          ///< @private
   void** prev;         ///< @private
   uint32_t seq;        ///< @private
   uint32_t typeId;     ///< @private
   uint32_t nor;        ///< Number of rows
   uint32_t noc;        ///< Number of columns
   int32_t* data;       ///< Marks (row by row)
} IntMatrix_t;

IntMatrix_t* imatAlloc(uint32_t nor, uint32_t noc);
int imatCompare(const IntMatrix_t* a, const IntMatrix_t* b);
IntMatrix_t* imatCreateFromBuffer(int32_t* buffer, uint32_t nor, uint32_t noc);
IntMatrix_t *imatDup(const IntMatrix_t* mat);
void imatFree(IntMatrix_t* mat);
IntMatrix_t* imatLoad(const char* fn);
IntMatrix_t* imatRead(MtxFile_t* file);
void imatSave(const IntMatrix_t* mat, const char* file_name);
void imatValidate(const struct MtxSourceLocation* sl, const IntMatrix_t* m);
void imatWrite(const IntMatrix_t* mat, MtxFile_t* f);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Polymorphic objects
////////////////////////////////////////////////////////////////////////////////////////////////////

void* objDup(void* a);
int objCanMultiply(void* a, void* b);
void objFree(void* a);
void* objInverse(void* a);
void* objLoad(const char* fn);
void objMul(void* a, void* b);
long objOrder(void* a);
void* objPower(void* a, int n);
int objSave(void* a, const char* fn);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Matrix representations
////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct {
   void* next;           ///< @private
   void** prev;          ///< @private
   uint32_t seq;         ///< @private
   uint32_t typeId;      ///< @private
   int NGen;
   Matrix_t** Gen;
} MatRep_t;

#define MR_COPY_GENERATORS  0x0001

void mrAddGenerator(MatRep_t* rep, Matrix_t* gen, int flags);
MatRep_t* mrAlloc(int ngen, Matrix_t* *gen, int flags);
int mrAreIsomorphic(const MatRep_t* repA, const MatRep_t* repB, Matrix_t* trans);
void mrChangeBasis(MatRep_t* rep, const Matrix_t* trans);
MatRep_t* mrChangeBasis2(const MatRep_t* rep, const Matrix_t* trans);
void mrValidate(const struct MtxSourceLocation* where, const MatRep_t* rep);
void mrFree(MatRep_t* rep);
MatRep_t* mrLoad(const char* basename, int ngen);
int mrSave(const MatRep_t* rep, const char* basename);
MatRep_t* mrTransposed(const MatRep_t* rep);

////////////////////////////////////////////////////////////////////////////////////////////////////
// The word generator
////////////////////////////////////////////////////////////////////////////////////////////////////

#define MTX_WG_MAXLEN 8

typedef struct {
   void* next;                  ///< @private
   void** prev;                 ///< @private
   uint32_t seq;                ///< @private
   uint32_t typeId;             ///< @private
   const MatRep_t* Rep;         ///< The representation
   Matrix_t* Basis[8];          ///< Products of the generators
   int N2[8];                   ///< Coefficients
   int* Description;            ///< Symbolic description of a word (binary)
   int buf[8][MTX_WG_MAXLEN + 1];  ///< internal
   int lastn2;
   char name[8 * (MTX_WG_MAXLEN + 1) + 1]; ///< Symbolic description of a word (text)
} WgData_t;

WgData_t* wgAlloc(const MatRep_t* rep);
int* wgDescribeWord(WgData_t* b, uint32_t n);
int wgFree(WgData_t* b);
Matrix_t* wgMakeWord(WgData_t* b, uint32_t n);
Matrix_t* wgMakeWord2(WgData_t* b, uint32_t n);
void wgMakeFingerPrint(WgData_t* b, uint32_t fp[6]);
const char* wgSymbolicName(WgData_t* b, long n);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Spin-up and split.
////////////////////////////////////////////////////////////////////////////////////////////////////

#define SF_SEED_MASK      0x0007  // Defines how to select seed vectors from the seed space
#define SF_FIRST          0x0001  // Use only the first basis vector
#define SF_EACH           0x0002  // Use each seed vector
#define SF_MAKE           0x0004  // Use all 1-dimensional subspaces of the seed space

#define SF_MODE_MASK      0x00F0  // Defines how to spin up and when to stop    
#define SF_SUB            0x0010  // Try finding a proper subspace
#define SF_CYCLIC         0x0020  // Try finding a cyclic vector (spins up to the whole space)
#define SF_COMBINE        0x0040  // Combine the spans

#define SF_STD            0x0100  // Spin up 'canonically' (standard basis)

#define SF_RESERVED_MASK  0xFFFFFE00

Matrix_t* quotientProjection(const Matrix_t* subspace, const Matrix_t* vectors);
Matrix_t* quotientAction(const Matrix_t* sub, const Matrix_t* gen);
Matrix_t* subspaceAction(const Matrix_t* sub, const Matrix_t* gen);

void split(const Matrix_t* subspace, const MatRep_t* rep, MatRep_t **sub, MatRep_t **quot);

int ConvertSpinUpScript(IntMatrix_t* script);


Matrix_t* spinup(const Matrix_t* seed, const MatRep_t* rep);
Matrix_t* spinupStandardBasis(
      IntMatrix_t** script, const Matrix_t* seed, const MatRep_t* rep, unsigned seedMode);
uint32_t spinupFindSubmodule(
      Matrix_t** basis,
      const Matrix_t* seed, const MatRep_t* rep, unsigned seedMode, uint32_t maxDimension);
uint32_t spinupFindCyclicVector(
      Matrix_t** basis,
      const Matrix_t* seed, const MatRep_t* rep, unsigned seedMode);

Matrix_t* spinupWithScript(const Matrix_t* seed, const MatRep_t* rep, const IntMatrix_t* script);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Seed vector generator
////////////////////////////////////////////////////////////////////////////////////////////////////

void svgMake(PTR vec, uint32_t number, const Matrix_t *basis);
int svgMakeNext(PTR vec, uint32_t* number, const Matrix_t* basis);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Miscellaneous algorithms
////////////////////////////////////////////////////////////////////////////////////////////////////

Matrix_t* matInsert_(Matrix_t* mat, const Poly_t* pol);
Matrix_t* matInsert(const Matrix_t* mat, const Poly_t* pol);
int IsSubspace(const Matrix_t* sub, const Matrix_t* space, int ngen);

Matrix_t* matTensor(const Matrix_t* m1, const Matrix_t* m2);
int MatrixToVector(const Matrix_t* mat, Matrix_t* vecs, int n);
Matrix_t* VectorToMatrix(Matrix_t* vecs, int n, int noc);
Matrix_t* TensorMap(Matrix_t* vec, const Matrix_t* a, const Matrix_t* b);

int StablePower(const Matrix_t* mat, int* pwr, Matrix_t **ker);
int StablePower_(Matrix_t* mat, int* pwr, Matrix_t **ker);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Polynomial factorization (Berlekamp algorithm)
////////////////////////////////////////////////////////////////////////////////////////////////////

FPoly_t* Factorization(const Poly_t* pol);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Characteristic and minimal polynomials
////////////////////////////////////////////////////////////////////////////////////////////////////

/// State for characteristic/minimal polynomial computation.
/// This structure is used by @ref charpolFactor / @ref minpolFactor.

struct CharpolState
{
   void* next;                  ///< @private
   void** prev;                 ///< @private
   uint32_t seq;                ///< @private
   uint32_t typeId;             ///< @private

   /// Defines which polynomial shall be calculated.
   enum CharpolMode {
      PM_CHARPOL,       ///< characteristic polynomial
      PM_MINPOL         ///< minimal polynomial
   } mode;

   long fl;      // field order
   long vsDim;   // vector space dimension
   long* piv;    // Pivot table
   char* ispiv;  // Pivot flags
   PTR mat;      // The matrix
   PTR A;        // Work space (for spin-up)
   PTR B;        // Work space II (coefficients)
   long dim;     // Dimension reached so far (sum of cyclic subspace dimensions)
   long n;       // Dimension of cyclic subspace
   long seed;    // Number of seed vector for the first cyclic subspace

   // Minimal polynomial on the current subspace dimension. Unused in PM_CHARPOL mode.
   Poly_t* partialMinPol;
};
typedef struct CharpolState Charpol_t;

Charpol_t* charpolStart(const Matrix_t* matrix, enum CharpolMode mode, long seed);
void charpolFree(struct CharpolState* state);

Poly_t* charpolFactor(Charpol_t* state);
FPoly_t* charpol(const Matrix_t* mat);

Poly_t* minpolFactor(Charpol_t* state);
FPoly_t* minpol(const Matrix_t* mat);
FPoly_t* minpolS(const Matrix_t* mat, long seed);

////////////////////////////////////////////////////////////////////////////////////////////////////
// GAP output support functions
////////////////////////////////////////////////////////////////////////////////////////////////////

const char* gapFelToString(FEL a);
const char* gapFelToString1(FEL a);
const char* gapFelToString2();
void gapFormatFel(StrBuffer_t* sb, FEL a);
void gapFormatPoly(StrBuffer_t* sb, const Poly_t* pol);
void gapFormatWord(StrBuffer_t* sb, const WgData_t* b, uint32_t n);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Submodule lattice functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#define MAXGEN 20       // Max. number of generators 
#define LAT_MAXCF 200   // Max. number of composition factors
#define MAXCYCL 30000   // Max. number of cyclic submodules
#define MAXDOTL 90000   // Max. number of dotted lines
#define MAXNSUB 20000   // Max. number of submodules

typedef struct {
   long dim;                    ///< Constituent dimension
   long num;                    ///< Constituent number (per dimension)
   int mult;                    ///< Multiplicity of the constituent
   char name[30];               ///< DimNum, e.g., "20b"
   uint32_t idWord;             ///< Identifying word number
   Poly_t* idPol;               ///< Identifying polynomial
   uint32_t peakWord;           ///< Peak word number
   Poly_t* peakPol;             ///< Peak polynomial
   long nmount;                 ///< Number of mountains
   long ndotl;                  ///< Number of dotted lines
   long spl;                    ///< Degree of splitting field
}
CfInfo;

typedef struct {
   void* next;                          ///< @private
   void** prev;                         ///< @private
   uint32_t seq;                        ///< @private
   uint32_t typeId;                     ///< @private
   char *baseName;                      ///< Module name
   int field;                           ///< Field order
   int NGen;                            ///< Number of generators
   int nCf;                             ///< Number of irred. constituents
   CfInfo Cf[LAT_MAXCF];                ///< Data for irred. constituents
   int NSocles;                         ///< Loewy length
   int* Socle;                          ///< Mult. of constituents in socles
   int NHeads;                          ///< Number of radical layers
   int* Head;                           ///< Mult. of constituents in Heads
} LatInfo_t;

LatInfo_t* latCreate(const char* baseName);
int latAddHead(LatInfo_t* li, int* mult);
int latAddSocle(LatInfo_t* li, int* mult);
const char* latCfName(const LatInfo_t* li, int cf);
void latDestroy(LatInfo_t* li);
LatInfo_t* latLoad(const char* basename);
void latSave(const LatInfo_t* li);

#define LAT_RG_INVERT           0x0001  // Invert generators
#define LAT_RG_TRANSPOSE        0x0002  // Transpose generators
#define LAT_RG_STD              0x0004  // Use standard form

MatRep_t* latReadCfGens(LatInfo_t* info, int cf, int flags);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Tensor condensation package
////////////////////////////////////////////////////////////////////////////////////////////////////

/// @addtogroup tp
/// @{

/// Tensor condensation state.
/// See also @ref tkReadInfo and @ref tkWriteInfo.

typedef struct {
   char *nameM;                 ///< Name of right factor
   char *nameN;                 ///< Name of left factor
   int dim;                     ///< Dimension of condensed module
   int nCf;                     ///< Number of relevant constituents
   int cfIndex[2][LAT_MAXCF];   ///< Constituent number in M/N
} TkData_t;

void tkReadInfo(TkData_t* tki, const char* name);
int tkWriteInfo(TkData_t* tki, const char* name);

/// @}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions operating on representations and vector spaces
////////////////////////////////////////////////////////////////////////////////////////////////////

int IsIsomorphic(const MatRep_t* rep1, const CfInfo* info1,
                 const MatRep_t* rep2, Matrix_t  **trans, int use_pw);
int makeEndomorphisms(const MatRep_t* rep, const Matrix_t* nsp,
                      Matrix_t* endo[]);
Matrix_t* HomogeneousPart(MatRep_t* m, MatRep_t* s, Matrix_t* npw,
                          const IntMatrix_t* op, int dimends);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Lattice drawing functions
////////////////////////////////////////////////////////////////////////////////////////////////////

// A lattice node.

typedef struct {
   double PosX, PosY;           // Position [0..1]
   unsigned long UserData;      // User-defined attributes
   int Layer;                   // Layer number
   double Score;                // Used in optimization
   int ScoreCount;
} LdNode_t;


// A lattice (nodes with x/y positions and parent/child relations).

typedef struct {
   int NNodes;
   LdNode_t* Nodes;
   int* IsSub;          // Incidence relation, <NNodes> * <NNodes> entries
   int* LayerNo;        // Layer numbers
   int NLayers;
} LdLattice_t;

#define LD_ISSUB(l,i,k) ((l)->IsSub[(i) * (l)->NNodes + (k)])

LdLattice_t* ldAlloc(int num_nodes);
int ldFree(LdLattice_t* l);
int ldAddIncidence(LdLattice_t* lat, int sub, int sup);
int ldSetPositions(LdLattice_t* l);

#endif  // !defined(_MEATAXE_H_)

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
