////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Initialization and clean up
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>
#include <limits.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @addtosection app
/// @{

static int isInitialized = 0;


////////////////////////////////////////////////////////////////////////////////////////////////////

static char libDir[250] = "";

int MtxOpt_UseOldWordGenerator = 0;



////////////////////////////////////////////////////////////////////////////////////////////////////

static void deriveDirectoryName(
      char *buf, size_t bufSize, const char *argv0, int strip, const char *suffix)
{
   if (argv0 != NULL && *argv0 == '/') {
      const char *end = argv0 + strlen(argv0);
      for (; strip > 0; --strip) {
         while (end > argv0 && end[-1] != '/')
            --end;
         if (end > argv0) --end;
      }
      if (strip == 0) {
         if (suffix != NULL && *suffix != 0)
             snprintf(buf, bufSize, "%.*s/%s", (int)(end - argv0), argv0, suffix);
         else
             snprintf(buf, bufSize, "%.*s", (int)(end - argv0), argv0);
         return;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void setDirectories(char* argv0)
{
   char* c;
   if ((c = getenv("MTXLIB")) != NULL) {
      snprintf(libDir, sizeof(libDir), "%s", c);
      return;
   } 
   if (libDir[0] == 0) {
      char* exePath = realpath(sysGetExecutableName(argv0), NULL);
      deriveDirectoryName(libDir, sizeof(libDir), exePath, 2, "lib");
      free(exePath);
   }
   if (libDir[0] == 0)
      snprintf(libDir, sizeof(libDir), ".");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int mtxIsBigEndian()
{
   static int value = -1;
   if (value == -1) {
      union {
         unsigned char c[sizeof(int)];
         unsigned u;
      } x;
      x.u = 1;
      value = (x.c[0] == 0);
   }
   return value;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the MeatAxe version.
/// This function can be called before @ref mtxInitLibrary.

const char* mtxVersion()
{
   static char version[200] = {0};
   if (version[0] == 0) {
      sysInit();
      char *v = version;
      char *vEnd = version + sizeof(version);
      v += snprintf(v, vEnd - v, "%s", MTX_VERSION);
      v += snprintf(v, vEnd - v, " I%u", (unsigned) sizeof(int) * 8);
      v += snprintf(v, vEnd - v, " L%u", (unsigned) sizeof(long) * 8);
      v += snprintf(v, vEnd - v, " %s", mtxIsBigEndian() ? "BE" : "LE");
      v += snprintf(v, vEnd - v, " ZZZ=%d ZZZVERSION=0x%x", MTX_ZZZ, MTX_ZZZVERSION);
   }
   return version;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// This function initializes the MeatAxe library including finite field arithmetic and file i/o
/// functions. It must be called before any other MeatAxe library function.
/// It is legal to call MtxInitLibrary() multiple times. Only the first call will actually do
/// anything.
/// An application that uses @ref MtxInitApplication need not call this function.
///
/// @a argv0 is the name of the process executable. It will be used to initialize directory names
/// such as @ref libDir, which have a default value relative to the executable directory. If the
/// program name is not known, the argument may be NULL or an empty string.

void mtxInitLibrary(char* argv0)
{
   if (isInitialized)
      return;
   isInitialized = 1;
   setDirectories(argv0);
   sysInit();
   if (sizeof(size_t) < sizeof(uint32_t))
      mtxAbort(MTX_HERE, "Unsupported platform");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the name of the MeatAxe library directory.
/// The returned name does not have a trailing slash, unless it is equal to "/".
/// This function fails and aborts the program if it is called before @ref mtxInitLibrary.

/// The library directory is determined as follows (in the order given here):
///
/// * If the "-L" option is used and the argument is not an empty string, the given directory
///   is added to the program's environment under the name MTXLIB.
/// * If the MTXLIB environment variable is defined and not an empty string, it is used as the
///   library directory.
/// * If MTXLIB is not defined or empty, the library directory is derived from the executable
///   directory by replacing the last path component with "lib". For example, if the program is
///   "/home/user1/mtx/bin/zcp", the derived library directory would be "/home/user1/mtx/lib".
/// * As a last resort, the current directory (".") is used.
///
/// There are no further checks whether the given directory exists and is usable.

const char *mtxLibraryDirectory()
{
   MTX_ASSERT(isInitialized);
   return libDir;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void mtxCleanupLibrary()
{
   memset(libDir, 0, sizeof(libDir));
   isInitialized = 0;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
