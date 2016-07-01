////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests extraction utility
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "check.h"

#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static test_Definition *tests = NULL;
static size_t nTests = 0;
static size_t maxTests = 0;
static const char *currentFileName;
static int currentLineNo;

////////////////////////////////////////////////////////////////////////////////////////////////////

MTX_PRINTF_ATTRIBUTE(1,2)
void fail(const char *message, ...)
{
   va_list args;
   va_start(args, message);
   fprintf(stderr, "ERROR: ");
   vfprintf(stderr, message, args);
   fprintf(stderr, "\n");
   va_end(args);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void printTests()
{
   printf("#include <check.h>\n");
   for (int i = 0; i < (int) nTests; ++i) {
      printf("extern test_F %s();\n", tests[i].name);
   }
   printf("test_Definition test_AllTests[] = {\n");
   for (int i = 0; i < (int) nTests; ++i) {
      printf("{%s,\"%s\", \"%s\", %d}%s\n",
	     tests[i].name,
	     tests[i].name,
	     tests[i].file,
	     tests[i].line,
	     i + 1 < (int) nTests ? "," : "");
   }
   printf("};\n");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void add(const char *name)
{
   if (maxTests >= nTests) {
      maxTests = 2 * nTests + 1;
      tests = (test_Definition*) realloc(tests, maxTests * sizeof(*tests));
   }
   test_Definition *t = tests + nTests;
   memset(t,0,sizeof(*t));
   t->name = name;
   t->file = currentFileName;
   t->line = currentLineNo;
   ++nTests;
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

static void addTest(const char *c)
{
   while (isspace(*c)) ++c;
   const char *son = c;
   while (isalnum(*c) || *c == '_') ++c;
   char *name = mkstr(son, c);

   add(name);
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
       if (!strncmp(lb, "test_F ", 7)) {
         addTest(lb + 7);
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
