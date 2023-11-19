////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Functions for Ephemeral Strings.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <string.h>


/// @addtogroup str
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////

#define MAX_EPHEMERAL_STRINGS 20

/// @private
struct TempBuffer {
   size_t index;
   char* string[MAX_EPHEMERAL_STRINGS];
};

#if defined(MTX_DEFAULT_THREADS)

static pthread_key_t tempKey;

static void createTempKey()
{
    pthread_key_create(&tempKey, NULL);
}

#else

static struct TempBuffer globalBuffer;

#endif


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Registers the given string for (eventual) deletion and returns the argument.
///
/// The argument must point to be a dynamically allocated block of memory. After strMakeEphemeral()
/// was called, the string is managed internally and must not be released or resized by the caller.
/// If the given string was already registerd for deletion, the function fails and aborts the
/// program.
///
/// The string buffer will be release automatically at a future call of strMakeEphemeral() in the
/// same thread. The implementation uses a FIFO buffer which keeps the last 20 strings alive.
/// Registering the 21st string releases and replaces the first string, and so on.
///
/// If multithreading is enabled, each thread uses an independent FIFO buffer. An ephemeral string
/// must only be used by the thread which created the string.

char* strMakeEphemeral(char* c)
{
   #if defined(MTX_DEFAULT_THREADS)
   static pthread_once_t createTmpKey = PTHREAD_ONCE_INIT;
   pthread_once(&createTmpKey, createTempKey);
   struct TempBuffer* temp = pthread_getspecific(tempKey);
   if (temp == NULL) {
      temp = ALLOC(struct TempBuffer);
      pthread_setspecific(tempKey, temp);
   }
   #else
   struct TempBuffer* temp = &globalBuffer;
   #endif

   for (size_t i = 0; i < MAX_EPHEMERAL_STRINGS; ++i) {
      if (temp->string[i] == c) {
         mtxAbort(MTX_HERE, "Multiple calls for the same string");
      }
   }

   if (temp->string[temp->index] != NULL) {
      memset(temp->string[temp->index], '!', strlen(temp->string[temp->index]));
      sysFree(temp->string[temp->index]);
   }
   temp->string[temp->index] = c;
   if (++temp->index == MAX_EPHEMERAL_STRINGS) {
       temp->index = 0;
   }
   return c;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
