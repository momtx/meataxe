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

////////////////////////////////////////////////////////////////////////////////////////////////////
/// This variable indicates if the MeatAxe library has been successfully
/// initialized (value 1) or not (value 0).

int Mtx_IsInitialized = 0;

int Mtx_IsBigEndian = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////

char MtxLibDir[250] = "";

int MtxOpt_UseOldWordGenerator = 0;

char MtxVersion[200] = "undefined";


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

   // Fallback to current directory.
   snprintf(buf, bufSize, ".");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void setDirectories(char* argv0)
{
   char* c;
   if ((c = getenv("MTXLIB")) != NULL) {
      snprintf(MtxLibDir, sizeof(MtxLibDir), "%s", c);
      return;
   } 

   char* path = realpath(argv0, NULL);
   if (MtxLibDir[0] == 0) {
      deriveDirectoryName(MtxLibDir, sizeof(MtxLibDir), path, 2, "lib");
   }
   free(path);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int MtxInitLibrary(char* argv0)
{
   if (Mtx_IsInitialized)
      return ZZZVERSION;
   Mtx_IsInitialized = 1;

   setDirectories(argv0);

   union {
      unsigned char c[sizeof(int)];
      unsigned u;
   } x;
   x.u = 1;
   Mtx_IsBigEndian = (x.c[0] == 0);

   sysInit();

   char *v = MtxVersion;
   char *vEnd = MtxVersion + sizeof(MtxVersion);
   v += snprintf(v, vEnd - v,"%s ", MTX_VERSION);
   if (sizeof(long) == 8) v += snprintf(v, vEnd - v, " L64");
   v += snprintf(v, vEnd - v, " %s", Mtx_IsBigEndian ? "BE" : "LE");
   v += snprintf(v, vEnd - v, " ZZZ=%d ZZZVERSION=0x%x", MTXZZZ, ZZZVERSION);



   return ZZZVERSION;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Terminate the library.
/// This function terminates the MeatAxe library. An application that uses
/// appFree() need not call this function.

void MtxCleanupLibrary()
{
}


void mtxSetLibraryDirectory(const char *dir)
{
   snprintf(MtxLibDir, sizeof(MtxLibDir), "%s", dir);
}


/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
