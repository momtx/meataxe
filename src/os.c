////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - OS dependent stuff
////////////////////////////////////////////////////////////////////////////////////////////////////

/// @defgroup os Operating System Interface
///
/// The MeatAxe is written for a UNIX-like operating environment and uses many functions of
/// the standard C library. To make the MeatAxe more portable between different operating
/// systems, some C library and system calls are accessed through wrapper functions. These
/// wrapper functions have names that begin with 'sys'. For example @c sysFree() is the wrapper
/// function for @c free().
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
// These sybols can be defined:
//
// OS_NO_CPU_TIME .......... no CPU times are available
// OS_TIMES_AND_SYSCONF	.... use times() and sysconf() instead of getrusage()
// OS_NO_ITIMER ............ no interval timers
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#if defined (_WIN32)

#pragma warning(disable: 4514)
#pragma warning(disable: 4201 4214 4115)
#include <windows.h>
#pragma warning(default: 4201 4214 4115)
#include <direct.h>

#else

#include <sys/stat.h>
#include <sys/time.h>
#include <sys/times.h>
#include <signal.h>
#include <unistd.h>
#include <sys/resource.h>

#endif

#include <errno.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

#if defined(OS_NO_CPU_TIME)
static time_t zinittime = 0;           /**< Start time of this process. */
#endif


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the smallest multiple of @a unit greater than or equal to @a x.

size_t sysPad(size_t x, size_t unit)
{
   size_t rem = x % unit;
   return rem == 0 ? x : x + (unit - rem);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// OS-specific initialization.
/// This function is called during library initialization. It performs any OS-specific
/// actions. Applications should never call this function directly. Use MtxInit() instead.

void sysInit()
{
#if defined(OS_NO_CPU_TIME)
   zinittime = time(NULL);
#endif
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// CPU time.
/// This function returns the CPU time used by the calling process in units of 1/10 seconds.
/// @see sysSetTimeLimit()
/// @return CPU time used.

long sysTimeUsed(void)
{
#if defined(_WIN32)

   FILETIME cr, ex, krnl, usr;
   LARGE_INTEGER k, u;
   GetProcessTimes(GetCurrentProcess(),&cr,&ex,&krnl,&usr);
   memcpy(&k,&krnl,sizeof(k));
   memcpy(&u,&usr,sizeof(u));
   return (long)((k.QuadPart + u.QuadPart) / 1000000);

#elif defined(OS_TIMES_AND_SYSCONF)

   struct tms t;
   static long clk_tck = 0;
   if (clk_tck == 0) {
      clk_tck = sysconf(_SC_CLK_TCK);
   }
   times(&t);
   return (long)((t.tms_utime + t.tms_stime) * 10 / clk_tck);

#elif defined(OS_NO_CPU_TIME)

   return (time(NULL) - zinittime) * 10;

#else

   static struct rusage ru;
   getrusage(RUSAGE_SELF,&ru);
   return ru.ru_utime.tv_sec * 10 + ru.ru_utime.tv_usec / 100000;

#endif

}


////////////////////////////////////////////////////////////////////////////////////////////////////

///
/// @fn sysSetTimeLimit(long)
/// Set CPU time limit.
/// This function sets a CPU time limit for the calling process. When the limit is exceeded
/// the process is killed.
/// @param nsecs CPU time limit in seconds.
/// @attention This function is not available on all platforms.

#ifdef _WIN32

DWORD CALLBACK Killer(void *x)
{
   DWORD ms = (long) x;
   Sleep(ms);
   mtxAbort(MTX_HERE,"%s",MTX_ERR_GAME_OVER);
   return 0;
}

void sysSetTimeLimit(long nsecs)
{
   DWORD id;
   CreateThread(NULL,0,Killer,(void *)(nsecs * 1000),0,&id);
}


#elif !defined(OS_NO_ITIMER)

void vtalarm(int i)
{
   mtxAbort(MTX_HERE,"%s",MTX_ERR_GAME_OVER);
}


void sysSetTimeLimit(long nsecs)
{
   struct itimerval tv;

   tv.it_interval.tv_sec = 0;
   tv.it_interval.tv_usec = 0;
   tv.it_value.tv_sec = nsecs;
   tv.it_value.tv_usec = 0;
   setitimer(ITIMER_VIRTUAL,&tv,NULL);
   signal(SIGVTALRM,vtalarm);
}


#else /* No interval timer, sorry... */

void sysSetTimeLimit(long nsecs)
{
}

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Opens a file. 
/// This function works like fopen() with the following differences:
/// * If the operation fails an error is raised, which normally aborts the program. If the
///   application has defined an error handler that doe snot abort, sysFopen() returns NULL on
///   error.
/// * The @a mode string can be extended by appending "::FLAGS", where FLAGS is
///   a colon-separated list of any of the following items:
///   * "lib" - Try to open the file in the library directory (see @ref mtxLibraryDirectory),
///     unless @a name starts with '/'. It this fails, try again using the file name as it is.
///     No errors are reported if the first attempt fails and the second attempt succeeds.
///   * "noerror" - Do not raise an error if the file cannot be opened, just return NULL.
///   For example: sysFopen("coeff7.txt", "r::lib:noerror")
///
/// @return A pointer to the open file or NULL on error.

FILE *sysFopen(const char *name, const char* mode)
{
   char sysMode[20];
   const char* mtxExt = strstr(mode, "::");
   int useLibDir = 0;
   int raiseError = 1;
   if (mtxExt != NULL) {
      snprintf(sysMode, sizeof(sysMode), "%.*s", (int)(mtxExt - mode), mode);
      const char* c = mtxExt + 2;
      while (1) {
         if (strncmp(c, "lib", 3) == 0) {
            useLibDir = 1;
            c += 3;
         }
         else if (strncmp(c, "noerror", 7) == 0) {
            raiseError = 0;
            c += 7;
         } else {
            break;
         }
         if (*c != ':') break;
         ++c;
      }
      if (*c != 0) {
         mtxAbort(MTX_HERE,"Invalid file mode %s", mode);
         return NULL;
      }
   } else {
      snprintf(sysMode, sizeof(sysMode), "%s", mode);
   }

   FILE *f = NULL;
   if (useLibDir && *name != '/') {
      const char* buf = strTprintf("%s/%s", mtxLibraryDirectory(), name);
      f = fopen(buf,sysMode);
   }
   if (f == NULL)
   {
      f = fopen(name,sysMode);
   }

   if (f == NULL && raiseError) {
      mtxAbort(MTX_HERE,"Cannot open %s (mode=%s): %s",name, sysMode, strerror(errno));
   }
   return f;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Set file pointer.
/// This function sets the file pointer to a given position. If pos is greater than or equal to
/// zero, it is interpreted as an absolute position (relative to start of file). If @c pos is
/// negative, the file pointer is moved to the end of file.
/// @see sysFseekRelative()
/// @param file File handle.
/// @param pos New position of file pointer.
/// @return 0 on success, nonzero otherwise.

int sysFseek(FILE *file, long pos)
{
   if (pos < 0) {
      return fseek(file,(long) 0,SEEK_END);
   } else {
      return fseek(file,pos,SEEK_SET);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Set file pointer relative to current position.
/// This function moves the file pointer by a given number of bytes, which may be positive
/// or negative.
/// @param file The file handle.
/// @param distance The number of bytes by which the file pointer shall be moved.
/// @return 0 on success, nonzero on error.
/// @see sysFseek()

int sysFseekRelative(FILE *file, long distance)
{
   return fseek(file,distance,SEEK_CUR);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Remove a file
/// This function deletes a file. On a UNIX system, sysRemoveFile() just calls remove().
/// If the file to be deleted does not exist or cannot be removed for some other reason,
/// run-time error error is generated.

int sysRemoveFile(const char *name)
{
   if (remove(name) != 0) {
      mtxAbort(MTX_HERE,"Cannot remove file '%s'",name);
      return -1;
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Removes a directory (which must be empty).
/// @see sysCreateDirectory().
/// @param name Name of the directory to be removed.
/// @return 0 on success, -1 on error.

int sysRemoveDirectory(const char *name)
{
   if (rmdir(name) != 0) {
      mtxAbort(MTX_HERE,"Cannot remove directory '%s'",name);
      return -1;
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Create a directory.
/// This function creates a new directory. If the directory cannot be created for some reason,
/// a run-time error is generated and the function returns -1.
/// @see sysRemoveDirectory()
/// @param name Name of the directory.
/// @return 0 on success, -1 on error.

int sysCreateDirectory(const char *name)
{
#ifdef _WIN32
   if (mkdir(name) != 0)
#else
   if (mkdir(name,0755) != 0)
#endif
   {
      mtxAbort(MTX_HERE,"Cannot create directory '%s'",name);
      return -1;
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Allocate memory.
/// This function works like @c malloc(), but the return value is never NULL, even if @a nbytes is
/// 0. The allocated memory region is initialized with zeroes.
/// @param nbytes Size of memory block to allocate.
/// @return Pointer to memory block.

void *sysMalloc(size_t nbytes)
{
   if (nbytes == 0) {
      nbytes = 1;       // avoid NULL return
   }
   void* x = calloc(1, nbytes);
   if (x == NULL) {
      mtxAbort(MTX_HERE,"Cannot allocate %lu bytes",(unsigned long) nbytes);
   }
   return x;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Resizes a memory block.
/// This function works like @c realloc() but handles zero-length blocks differently (namely, by
/// allocating 1 byte instead) to avoid problems with broken @c realloc() implementations.
/// Note: if a buffer is enlarged, the added memory region is not initialized with zeroes!
/// @param buf Pointer to the memory block.
/// @param nbytes Desired new size.
/// @return  Pointer to resized memory block or NULL on error.

void *sysRealloc(void *buf, size_t nbytes)
{
   void *x;
   if (nbytes == 0) {
      nbytes = 1;
   }
   x = realloc(buf,nbytes);
   if (x == NULL) {
      mtxAbort(MTX_HERE,"Cannot reallocate %lu bytes: %s",(long) nbytes, strerror(errno));
   }
   //printf("realloc(0x%p, %lu)=0x%p\n", buf, (unsigned long) nbytes, x);
   return x;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// This function works like @c free().

void sysFree(void *x)
{
   //printf("free(0x%p)\n", x);
   if (x != NULL) {
      free(x);
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Get process id.
/// This function returns a number which uniquely identifies the calling process on the local
/// system. The exact meaning of this number depends on the operating system. In an UNIX
/// environment, it is the process id (PID).
/// @return Process id.

int sysGetPid()
{
   int pid;
#ifdef _WIN32
   pid = (int) GetCurrentProcessId();
#else
   pid = (int) getpid();
#endif
   return pid;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

const char* sysGetExecutableName(const char* argv0)
{
   static char *result = NULL;
   if (result != NULL)
      return result;

   if (strchr(argv0,'/') != NULL)
      return argv0;
   const char* p = getenv("PATH");

   while (1) {
      while (*p == ':') ++p;
      if (*p == 0) break;
      const char *bop = p;
      while (*p != 0 && *p != ':') ++p;
      int len = p - bop;

      const char* exeName = strTprintf("%.*s/%s", len, bop, argv0);
      if (access(exeName, X_OK) == 0) {
         result = strdup(exeName);
         return result;
      }
   }
   return argv0;
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
