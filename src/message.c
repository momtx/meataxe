////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Messages.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>

/// @addtogroup app
/// @{

/// @def MTX_XLOGE(sb)
/// Create a complex log message (level MTX_LOG_ERROR).
/// The MTX_XLOG? macros can be used to write complex log messages with minimal overhead if the
/// message is disabled by the current log threshold. The macro must be followed by a block
/// statement which has access to a string buffer (see @ref StrBuffer) under the name @a sb.
/// Inside the block statement you can use @ref sbAppend, @ref sbPrintf or any other formatting
/// function that works with a string buffer target. The block will only be executed if the message
/// is permitted by the current log threshold.
///
/// Example:
/// <pre>
/// MTX_XLOG_E(message) {
///    sbAppend(message, "polynomial=");
///    Poly_t* polynomial = costlyFunctionCalculatingThePolynomial();
///    polFormat(message, polynomial);
///    polFree(polynomial);
///    sbPrintf(message, ", multiplicity=%d", multiplicity);
/// }
/// </pre>
/// 
/// As you can see in the example, the message buffer is managed automatically and must not be
/// released. Any other dynamically allocated objects must be released inside the block.

/// @def MTX_XLOGW(sb)
/// Create a complex log message (level MTX_LOG_WARNING). See @ref MTX_XLOGE.

/// @def MTX_XLOGI(sb)
/// Create a complex log message (level MTX_LOG_INFO). See @ref MTX_XLOGE.

/// @def MTX_XLOGD(sb)
/// Create a complex log message (level MTX_LOG_DEBUG). See @ref MTX_XLOGE.

/// @def MTX_XLOG2(sb)
/// Create a complex log message (level MTX_LOG_DEBUG2). See @ref MTX_XLOGE.

static int defaultThreshold = 0;
static int logThreshold = 0;
static int fmtHash = 0;
static int fmtLevel = 0;
static int fmtThread = 0;
static int fmtTime = 0;

static FILE* logFile = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////

void logSetDefaultThreshold(int level)
{
   if (level > MTX_LOG_DEBUG2)
      defaultThreshold = MTX_LOG_DEBUG2;
   else if (level < MTX_LOG_ERROR)
      defaultThreshold = MTX_LOG_ERROR;
   else
      defaultThreshold = level;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int logGetDefaultThreshold()
{
   return defaultThreshold;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int init(int level)
{
   if (logFile == NULL) {
      logFile = stdout;
      logThreshold = defaultThreshold;
      fmtLevel = 0;
      fmtHash = 0;
      fmtTime = 0;
      fmtThread = 0;
   }
   if (level > logThreshold)
      return -1;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void startLine(StrBuffer* sb, int level)
{
   if (fmtHash)
      sbAppend(sb, "#");
   if (fmtTime) {
      struct timeval now;
      gettimeofday(&now, NULL);
      struct tm tm_now;
      localtime_r(&now.tv_sec, &tm_now);
      sbPrintf(sb, "%s%04d-%02d-%02d %02d:%02d:%02d.%03d ",
         sb->size == 0 ? "" : " ",
         tm_now.tm_year + 1900,
         tm_now.tm_mon + 1,
         tm_now.tm_mday,
         tm_now.tm_hour,
         tm_now.tm_min,
         tm_now.tm_sec,
         (int)now.tv_usec / 1000);
   }
   if (fmtLevel) {
      const char* lvl = "d ";
      switch (level) {
         case MTX_LOG_ERROR:   lvl = "E "; break;
         case MTX_LOG_WARNING: lvl = "W "; break;
         case MTX_LOG_INFO:    lvl = "I "; break;
         case MTX_LOG_DEBUG:   lvl = "D "; break;
         default:
            lvl = (level < MTX_LOG_ERROR) ? "E " : "d ";
            break;
      }
      sbAppend(sb, lvl);
   }
   if (fmtThread) {
      const char* thread = pexLogPrefix();
      if (thread)
         sbAppend(sb, thread);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

StrBuffer* logStart(int level)
{
   if (init(level) != 0)
      return NULL;
   StrBuffer* sb = sbAlloc(100);
   startLine(sb, level);
   return sb;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns 1 if the given log level is enabled and 0 otherwise.

int logEnabled(int level)
{
   return level <= logThreshold;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void logBuffered(StrBuffer* buf)
{
   sbAppend(buf, "\n");
   fwrite(buf->data, 1, buf->size, logFile);
   sbFree(buf);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void logPrintf(int level, const char* msg, ...)
{
   StrBuffer* buffer = logStart(level);
   if (buffer == NULL)
      return;
   va_list args;
   va_start(args,msg);
   sbVprintf(buffer, msg, args);
   va_end(args);
   logBuffered(buffer);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void setFile(const char* begin, const char* end)
{
   if (logFile != NULL) {
      fclose(logFile);
      logFile = NULL;
   }

   const char *mode = "w";
   if (*begin == '+') {
      mode = "a";
      ++begin;
   }
   char* fileName = strRange(begin, end);
   if (strcmp(fileName, "stdout") == 0 || *fileName == 0) {
      logFile = stdout;
   }
   else if (strcmp(fileName, "stderr") == 0) {
      logFile = stderr;
   }
   else {
      logFile = sysFopen(fileName, mode);
   }
   sysFree(fileName);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int setThreshold(const char* begin, const char *end)
{
   if (begin == end) {
      logThreshold = defaultThreshold;
   }
   else if (strCompareRange(begin, end, "info", NULL) == 0) {
      logThreshold = MTX_LOG_INFO;
   }
   else if (strCompareRange(begin, end, "error", NULL) == 0) {
      logThreshold = MTX_LOG_ERROR;
   }
   else if (strCompareRange(begin, end, "warning", NULL) == 0) {
      logThreshold = MTX_LOG_WARNING;
   }
   else if (strCompareRange(begin, end, "debug", NULL) == 0) {
      logThreshold = MTX_LOG_DEBUG;
   }
   else if (strCompareRange(begin, end, "debug2", NULL) == 0) {
      logThreshold = MTX_LOG_DEBUG2;
   }
   else
      return -1;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int setFormat(const char* spec)
{
   if (*spec == 0 || strcmp(spec, "default") == 0) {
      fmtLevel = 1;
      fmtTime = 0;
   }
   else if (strcmp(spec, "none") == 0) {
      fmtTime = 0;
      fmtThread = 0;
      fmtLevel = 0;
   }
   else if (strcmp(spec, "full") == 0) {
      fmtLevel = 1;
      fmtThread = 1;
      fmtTime = 1;
   } else {
      return -1;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Configures log output.
/// @a spec format: «file»:«threshold»[:«format»]
///
/// If logInit() is not called, logging works with a default configuration of
/// "stdout:debug:default".

void logInit(const char* spec)
{
   const char* const sep1 = strchr(spec,':');
   const char* const sep2 = sep1 ? strchr(sep1 + 1, ':') : NULL;
   if (sep1 == NULL)
      mtxAbort(MTX_HERE, "Invalid log specification (missing ':'): \"%s\"", spec);
   setFile(spec, sep1);
   if (setThreshold(sep1 + 1, sep2) != 0) {
      mtxAbort(MTX_HERE, "Invalid log specification (unknown level): \"%s\"", spec);
   }
   if (setFormat(sep2 ? sep2 + 1 : "default") != 0) {
      mtxAbort(MTX_HERE, "Invalid log specification (unknown format): \"%s\"", spec);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void logPrepareForAbort()
{
   if (logFile == NULL)
      logFile = stderr;
   logThreshold = MTX_LOG_INFO;
   fmtThread = 1;
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
