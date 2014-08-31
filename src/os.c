/* ============================= C MeatAxe ==================================
   File:        $Id: os.c,v 1.4 2007-09-09 21:38:11 mringe Exp $
   Comment:     OS dependent stuff
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


/** @defgroup os Operating System Interface
  * @{
  * @details
  * The MeatAxe is written for a UNIX-like operating environment and uses many functions of
  * the standard C library. To make the MeatAxe more portable between different operating
  * systems, some C library and system calls are accessed through wrapper functions. These
  * wrapper functions have names that begin with 'Sys'. For example @c SysFree() is the wrapper
  * function for @c free().
  **/


/* --------------------------------------------------------------------------
   These sybols can be defined:

   OS_NO_CPU_TIME .......... no CPU times are available
   OS_TIMES_AND_SYSCONF	.... use times() and sysconf() instead of getrusage()
   OS_NO_ITIMER ............ no interval timers
   -------------------------------------------------------------------------- */


/* --------------------------------------------------------------------------
   Include files
   -------------------------------------------------------------------------- */

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


#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdarg.h>


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

#if defined(OS_NO_CPU_TIME)
time_t zinittime = 0;		/**< Start time of this process. */
#endif


MTX_DEFINE_FILE_INFO

/* ------------------------------------------------------------------
   fopen() modes
   ------------------------------------------------------------------ */

static char *fmodes[7] = { NULL,"rt","wt","at","rb","wb","ab" };



/**
 ** OS-specific initialization.
 ** This function is called during library initialization. It performs any OS-specific
 ** actions. Applications should never call this function directly. Use MtxInit() instead.
 **/

void SysInit()
{
#if defined(OS_NO_CPU_TIME)
    zinittime = time(NULL);
#endif
}



/**
 ** CPU time.
 ** This function returns the CPU time used by the calling process in units of 1/10 seconds.
 ** @see SysSetTimeLimit()
 ** @return CPU time used.
 **/

long SysTimeUsed(void)
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
    if (clk_tck == 0) 
	clk_tck = sysconf(_SC_CLK_TCK);
    times(&t);
    return ((long)((t.tms_utime + t.tms_stime) * 10 / clk_tck ));

#elif defined(OS_NO_CPU_TIME)

    return (time(NULL) - zinittime) * 10;

#else

    static struct rusage ru;
    getrusage(RUSAGE_SELF,&ru);
    return (ru.ru_utime.tv_sec * 10 + ru.ru_utime.tv_usec / 100000);

#endif


}



/**
 ** @fn SysSetTimeLimit(long)
 ** Set CPU time limit.
 ** This function sets a CPU time limit for the calling process. When the limit is exceeded
 ** the process is killed.
 ** @param nsecs CPU time limit in seconds.
 ** @attention This function is not available on all platforms.
 **/

#ifdef _WIN32

DWORD CALLBACK Killer(void *x)
{
    DWORD ms = (long) x;
    Sleep(ms);
    MTX_ERROR1("%E",MTX_ERR_GAME_OVER);
    return 0;
}

void SysSetTimeLimit(long nsecs)
{
    DWORD id;
    CreateThread(NULL,0,Killer,(void *)(nsecs * 1000),0,&id);
}

#elif !defined(OS_NO_ITIMER)

void vtalarm(int i)

{
    MTX_ERROR1("%E",MTX_ERR_GAME_OVER);
}

void SysSetTimeLimit(long nsecs)
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

void SysSetTimeLimit(long nsecs)

{
}

#endif



/**
 ** Open a file.
 ** This function opens a file, like @c fopen(). The second argument, must be one of the
 ** predfined constants FM_READ (open for reading), FM_CREAT (create a new file and open for
 ** writing, or FM_APPEND (append to existing file or create a new file).
 ** Additional flags may be or'ed to the mode:
 ** @par FM_LIB
 ** If the file does not exist in the current directory, look in the library directory.
 ** The library directory is defined either by the environment variable @c MTXLIB, or at
 ** compile-time by the macro @c MTXLIB.
 ** @par FM_TEXT
 ** Open in text mode. This flag must be used on some systems (e.g., MS-DOS) to open text files.
 ** By default, files are assumed to contain binary data.
 ** @par FM_NOERROR
 ** Do not generate an error if the file does not exist.
 ** @see FfReadHeader() FfWriteHeader()
 ** @return A pointer to the open file or NULL on error.
 **/

FILE *SysFopen(const char *name, int mode)
{
    char buf[300];
    int m;
    FILE *f;

    m = mode & 0x0F;			/* Append, read or create? */
    if ((mode & FM_TEXT) == 0) m += 3;	/* Binary mode */
    if (m < 1 || m > 6 || (mode & 0x0F) == 0)
    {
	MTX_ERROR1("Invalid file mode %d",mode);
	return NULL;
    }
    f = fopen(name,fmodes[m]);
    if (f != NULL) 
	return f;

    /* Search library directory
       ------------------------ */
    if ((mode & FM_LIB) != 0) 
    {
	strcpy(buf,MtxLibDir);
	strcat(buf,"/");
	strcat(buf,name);
	f = fopen(buf,fmodes[m]);
    }

    /* Error handling
       -------------- */
    if (f == NULL && (mode & FM_NOERROR) == 0)
	MTX_ERROR1("%s: %S",name);

    return f;
}


/** Set file pointer.
 ** This function sets the file pointer to a given position. If pos is greater than or equal to
 ** zero, it is interpreted as an absolute position (relative to start of file). If @c pos is
 ** negative, the file pointer is moved to the end of file.
 ** @see SysFseekRelative(), FfSeekRow()
 ** @param file File handle.
 ** @param pos New position of file pointer.
 ** @return 0 on success, nonzero otherwise.
 **/

int SysFseek(FILE *file, long pos)
{
    if (pos < 0)
	return fseek(file,(long) 0,SEEK_END);
    else
	return fseek(file,pos,SEEK_SET);
}



/**
 ** Remove a file
 ** This function deletes a file. On a UNIX system, SysRemoveFile() just calls remove().
 ** If the file to be deleted does not exist or cannot be removed for some other reason,
 ** run-time error error is generated.
 **/ 

int SysRemoveFile(const char *name)
{
    if (remove(name) != 0)
    {
	MTX_ERROR1("Cannot remove file '%s'",name);
	return -1;
    }
    return 0;
}




/**
 ** Remove a directory.
 ** This function removes the specified directory.
 ** @see SysCreateDirectory().
 ** @param name Name of the directory.
 ** @return 0 on success, -1 on error.
 **/

int SysRemoveDirectory(const char *name)
{
    if (rmdir(name) != 0)
    {
	MTX_ERROR1("Cannot remove directory '%s'",name);
	return -1;
    }
    return 0;
}





/**
 ** Create a directory.
 ** This function creates a new directory. If the directory cannot be created for some reason,
 ** a run-time error is generated and the function returns -1.
 ** @see SysRemoveDirectory()
 ** @param name Name of the directory.
 ** @return 0 on success, -1 on error.
 **/

int SysCreateDirectory(const char *name)
{
#ifdef _WIN32
    if (mkdir(name) != 0)
#else
    if (mkdir(name,0755) != 0)
#endif
    {
	MTX_ERROR1("Cannot create directory '%s'",name);
	return -1;
    }
    return 0;
}


/**
 ** Set file pointer relative to current position.
 ** This function moves the file pointer by a given number of bytes, which may be positive
 ** or negative.
 ** @param file The file handle.
 ** @param distance The number of bytes by which the file pointer shall be moved.
 ** @return 0 on success, nonzero on error.
 ** @see SysFseek(), FfSeekRow()
 **/

int SysFseekRelative(FILE *file, long distance)
{
    return fseek(file,distance,SEEK_CUR);
}


/**
 ** Allocate memory.
 ** This function works like @c malloc(), but the return value is never 0, even when the function
 ** was called with a 0 argument.
 ** @param nbytes Size of memory block to allocate.
 ** @return Pointer to memory block or NULL on error.
 **/

void *SysMalloc(size_t nbytes)
{
    void *x;

    if (nbytes == 0) 
	nbytes = 1;
    x = malloc(nbytes);
    if (x == NULL)
	MTX_ERROR1("Cannot allocate %l bytes: %S",(long int) nbytes);
    return x;
}


/**
 ** Resize memory block.
 ** This function works like @c realloc() but handles zero-length blocks differently (namely, by
 ** allocating 1 byte instead) to avoid problems with broken @c realloc() implementations.
 ** @param buf Pointer to memory block.
 ** @param nbytes Desired new size.
 ** @return  Pointer to resized memory block or NULL on error.
 **/

void *SysRealloc(void *buf, size_t nbytes)
{
    void *x;
    if (nbytes == 0) 
	nbytes = 1;
    x = realloc(buf,nbytes);
    if (x == NULL)
	MTX_ERROR1("Cannot reallocate %l bytes: %S",(long) nbytes);
    return x;
}


/**
 ** Free memory block.
 ** This function works like @c free() but checks if the argument is not NULL. Otherwise, an 
 ** appropriate error message is generated.
 ** @param x Pointer to the memory block.
 **/

void SysFree(void *x)
{
    if (x == NULL)
	MTX_ERROR("Attempt to free() NULL pointer");
    else
	free(x);
}


/** Get process id.
 ** This function returns a number which uniquely identifies the calling process on the local
 ** system. The exact meaning of this number depends on the operating system. In an UNIX
 ** environment, it is the process id (PID).
 ** @return Process id.
 **/

int SysGetPid()
{
    int pid;
#ifdef _WIN32
    pid = (int) GetCurrentProcessId();
#else
    pid = (int) getpid();
#endif
    return pid;
}


/** @}
 **/

