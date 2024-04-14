////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests extraction utility
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "testing.h"

#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char *currentFileName;
static int currentLineNo;

/// @private
struct FoundTest {
    struct FoundTest* prev;
    char* name;
    int flags;
};

struct FoundTest* foundTests = NULL;



////////////////////////////////////////////////////////////////////////////////////////////////////

MTX_PRINTF(1,2)
void fail(const char *message, ...)
{
   va_list args;
   va_start(args, message);
   if (currentFileName != NULL) {
       fprintf(stderr, "%s:%d:: error:", currentFileName, currentLineNo);
   }
   vfprintf(stderr, message, args);
   fprintf(stderr, "\n");
   va_end(args);
   exit(1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int strSkipExact(const char**rpp, const char* pattern)
{
   const char *rp = *rpp;
   while (*pattern != 0) {
      if (*rp == *pattern) {
         ++rp;
      } else {
         return 0;
      }
      ++pattern;
   }
   *rpp = rp;
   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int strSkip(const char**rpp, const char* pattern)
{
   const char *rp = *rpp;
   while (*pattern != 0) {
      if (*pattern == ' ') {
         // Zero or more spaces
         while (isspace(*rp)) ++rp;
      }
      else if (*rp == *pattern) {
         ++rp;
      }
      else
         return 0;
      ++pattern;
   }
   *rpp = rp;
   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Returns a copy of «c» up to (not including) «cEnd».
// If «cEnd» is NULL, the full string is copied.

static char *strCopyRange(const char* c, const char* cEnd)
{
   const size_t len = cEnd ? cEnd - c : strlen(c);
   char *s = (char*) malloc(len + 1);
   memcpy(s, c, len);
   s[len] = 0;
   return s;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Returns a copy of the identifier at the read pointer position or NULL on error.

static char* strParseIdentifier(const char**rpp)
{
    const char *rp = *rpp;

    if (!isalpha(*rp)) return NULL;
    const char *start = rp;
    while (isalnum(*rp) || *rp == '_') ++rp;
    char *id = strCopyRange(start, rp);
    *rpp = rp;
    return id;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void printTests()
{
   printf("#include <testing.h>\n");
   for (const struct FoundTest* t = foundTests; t; t = t->prev)
   {
      printf("extern TstResult %s(", t->name);
      if (t->flags & TST_FLAG_PER_FIELD)
         printf("int q");
      printf(");\n");
   }

   printf("struct TstFoundTest foundTests[] = {\n");
   for (const struct FoundTest* t = foundTests; t; t = t->prev)
   {
      const char* localSuffix = strstr(t->name, "__");
      char* name = localSuffix ? strCopyRange(t->name, localSuffix) : NULL;
      printf("{%s, 0x%x, \"%s\"},\n", t->name, t->flags, name ? name : t->name);
      free(name);
   }
   printf("{NULL,0}\n");
   printf("};\n");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static char *mkstr(const char *begin, const char *end)
{
   const size_t len = end ? end - begin : strlen(begin);
   char *s = (char*) malloc(len + 1);
   memcpy(s,begin,len);
   s[len] = 0;
   return s;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

struct FoundTest* addFoundTest(char* name, unsigned int flags)
{
    struct FoundTest* t = (struct FoundTest*) calloc(1, sizeof(struct FoundTest));
    t->name = name;
    t->flags = flags;
    t->prev = foundTests;
    foundTests = t;
    return t;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void registerTest(const char *rp)
{
    char *testName = strParseIdentifier(&rp);
    if (testName == NULL) fail("Missing test name");
    if (!strSkip(&rp, " (")) fail("Missing \"(\" after test name");
    int flags = 0;
    if (strSkip(&rp,"int q )"))
	flags |= TST_FLAG_PER_FIELD;
    else if (!strSkip(&rp, " )"))
	fail("Missing \")\" in test function");
    addFoundTest(testName, flags);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void extract(const char *fileName)
{
    currentFileName = mkstr(fileName, 0);
    FILE *f = fopen(fileName,"r");
    if (f == NULL) {
	fail("Cannot open \"%s\": %s", fileName, strerror(errno));
    }
    char lb[10000];
    currentLineNo = 0;
    while (fgets(lb,sizeof(lb),f) != 0) {
	++currentLineNo;
	const char* rp = lb;
	if (strSkipExact(&rp, "TstResult ")) {
	    registerTest(rp);
	}
    }
    fclose(f);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
   for (int i = 1; i < argc; ++i) {
       extract(argv[i]);
   }
   printTests();
   return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
