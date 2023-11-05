////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Convert a matrix or permutation from ASCII (readable)
////////////////////////////////////////////////////////////////////////////////////////////////////

#define MAXLINE 4000    // max. input line size

#include "meataxe.h"
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static int GrpLibFormat = 0;    /* File is in Group Library Format */
static FILE *src = NULL;                // Input file
static FILE *out;                       // Output file
static int inpLineNo = 0;               // Input line number
static char lbuf[MAXLINE] = {0};        // Input line buffer
static char *lptr = lbuf;               // Read pointer
static const char *inpname = "[stdin]";
static const char *outname = "";
static int MemberCount = 0;             // Number of members

static MtxApplicationInfo_t AppInfo = {
   "zcv","Convert Text to Binary Format",
   "SYNTAX\n"
   "    zcv <Inp> <Out>\n"
   "\n"
   "ARGUMENTS\n"
   "    <Inp> ................... Input file (MeatAxe text format). '-' for stdin\n"
   "    <Out> ................... Output file (MeatAxe binary format)\n"
   "\n"
   "FILES\n"
   "    <Inp> ................... I Text file\n"
   "    <Out> ................... O Binary file\n"
};

static MtxApplication_t *App = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read next input line, skip empty lines and strip comments.
/// Returns 1 on success, 0 on end-of-file.

static int tryReadLine()
{
   while (1) {
      lbuf[0] = 0;
      if (fgets(lbuf,sizeof(lbuf),src) == NULL) {
         if (feof(src)) 
            return 0;
         mtxAbort(MTX_HERE,"Error reading file");
      }
      ++inpLineNo;

      // Skip comment lines (starting with '#')
      if (lbuf[0] == '#')
         continue;

      // Trim leading and trailing white space, skip empty lines.
      lptr = lbuf;
      while (isspace(*lptr))
         ++lptr;
      if (*lptr == 0)
         continue;
      char *c = lptr + strlen(lptr);
      if (c == lbuf + sizeof(lbuf) && c[-1] != '\n')
          mtxAbort(MTX_HERE,"Line too long");
      while (c > lptr && isspace(c[-1]))
         --c;
      *c = 0;
      return 1;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Like tryReadLine, but fails on end-of-file.

static void readLine()
{
   if (!tryReadLine())
      mtxAbort(MTX_HERE,"Unexpected end of file");
}

const char* provideInputFilePosition(void* userData)
{
   static char buffer[300];
   snprintf(buffer, sizeof(buffer), "Reading %s, line %d, column %d",
         inpname, inpLineNo, (int)(lptr - lbuf) + 1);
   return buffer;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

static void assertEndOfLine()
{
   while (isspace(*lptr))
      ++lptr;
   if (*lptr != 0)
      mtxAbort(MTX_HERE, "Unexpected trailing characters");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

/// Read field element

static FEL readFel()
{
   // Skip to beginning of next number.
   while (1) {
      while (isspace(*lptr)) ++lptr;
      if (*lptr != 0) 
         break;
      readLine();
   }

   // Parse number
   char *rp = lptr;
   uint32_t val = 0;
   if (!isdigit(*lptr))
      mtxAbort(MTX_HERE, "Bad input: expected digit, found 0x%02x", (unsigned char) *lptr);
   while (isdigit(*rp) && val < 0xFFFF) {
      val = val * 10 + (*rp - '0');
      ++rp;
      // Allow packed format without spaces for GF(q), q < 10
      if (ffOrder < 9)
         break;
   }
   if (val >= ffOrder) {
      mtxAbort(MTX_HERE, "Bad input: %lu not in GF(%lu)",
            (unsigned long) val, (unsigned long) ffOrder);
   }
   lptr = rp;
   return ffFromInt(val);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////

/// Read integer number

static int UNSIGNED = 1;
static int SIGNED = 2;

static uint32_t read32(const int type)
{
   // Skip to beginning of next number.
   while (1) {
      while (isspace(*lptr)) ++lptr;
      if (*lptr != 0) 
         break;
      readLine();
   }

   // Parse number
   int minus = 0;
   if (type == SIGNED) {
      if (*lptr == '-') {
         minus = 1;
         ++lptr;
      }
   }
   char *rp = lptr;
   uint64_t val = 0;
   const uint64_t maxVal = (type == SIGNED) ? (minus ? 0x80000000 : 0x7FFFFFFF) : 0xFFFFFFFF;
   if (!isdigit(*rp))
      mtxAbort(MTX_HERE, "Bad input: expected digit, found 0x%02x", (unsigned char) *rp);
   while (isdigit(*rp)) {
      val = val * 10 + (*rp - '0');
      if (val > maxVal)
          mtxAbort(MTX_HERE, "Number out of range [-2^32,2^32-1]");
      ++rp;
   }
   if (*rp != 0 && !isspace(*rp))
      mtxAbort(MTX_HERE, "Malformed number");

   lptr = rp;
   return (uint32_t)(minus ? -val : val);
}

static int32_t parse32s()
{
   return (int32_t) read32(SIGNED);
}

static uint32_t parse32u()
{
   return read32(UNSIGNED);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void WriteHeader(uint32_t a, uint32_t b, uint32_t c)
{
   uint32_t hdr[3] = {a, b, c};
   sysWrite32(out,hdr,3);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

//static const char* strPrefix(const char* s, const char* prefix)
//{
//   if (s == NULL || prefix == NULL)
//      return NULL;
//   while (*prefix != 0 && *s == *prefix) {
//      ++s;
//      ++prefix;
//   }
//   return (*prefix == 0) ? s : NULL;
//}
//
////////////////////////////////////////////////////////////////////////////////////////////////////

int tryParseLiteral(const char *text)
{
   size_t textLen = strlen(text);
   if (strncmp(lptr, text, textLen) != 0)
      return 0;
   char *end = lptr + textLen;
   if (*end != 0 && !isspace(*end))
      return 0;
   lptr = end;
   return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int tryParseHeaderField(uint32_t* var, const char* prefix)
{
   while (isspace(*lptr)) ++lptr;
   const size_t plen = strlen(prefix);
   if (strncmp(lptr, prefix, plen) != 0)
      return 0;
   lptr += plen;
   char* end = NULL;
   unsigned long long value = strtoull(lptr, &end, 10);
   if (*lptr != 0 && end != lptr && value <= 0xFFFFFFFFU) {
      *var = (uint32_t) value;
      lptr = end;
      return 1;
   }
   mtxAbort(MTX_HERE, "Invalid number after \"%s\"", prefix);
   return 0;
}
         
////////////////////////////////////////////////////////////////////////////////////////////////////

static int tryParseHeaderField32s(int32_t* var, const char* prefix)
{
   while (isspace(*lptr)) ++lptr;
   const size_t plen = strlen(prefix);
   if (strncmp(lptr, prefix, plen) != 0)
      return 0;
   lptr += plen;
   char* end = NULL;
   long long value = strtoll(lptr, &end, 10);
   if (*lptr != 0 && end != lptr && value <= 2147483647LL && value >= -2147483648LL) {
      *var = (int32_t) value;
      lptr = end;
      return 1;
   }
   mtxAbort(MTX_HERE, "Invalid number after \"%s\"", prefix);
   return 0;
}
         
////////////////////////////////////////////////////////////////////////////////////////////////////

static void convertMatrix()
{
   uint32_t fl = MTX_NVAL;
   uint32_t nor = MTX_NVAL;
   uint32_t noc = MTX_NVAL;
   while (1) {
      if (tryParseHeaderField(&fl, "field=")) {continue;}
      if (tryParseHeaderField(&nor, "nor=")) {continue;}
      if (tryParseHeaderField(&nor, "rows=")) {continue;}
      if (tryParseHeaderField(&noc, "noc=")) {continue;}
      if (tryParseHeaderField(&noc, "cols=")) {continue;}
      if (*lptr == 0) {break;}
      mtxAbort(MTX_HERE, "Unknown header field");
   }
   if (nor == MTX_NVAL) {mtxAbort(MTX_HERE, "Missing header field \"rows\".");}
   if (noc == MTX_NVAL) {mtxAbort(MTX_HERE, "Missing header field \"cols\".");}
   if (fl == MTX_NVAL) {mtxAbort(MTX_HERE, "Missing header field \"field\".");}

   MESSAGE(0, ("%dx%d matrix over GF(%d)\n", nor, noc, fl));
   ffSetField(fl);
   PTR m1 = ffAlloc(1, noc);
   WriteHeader(fl, nor, noc);

   for (uint32_t i = 0; i < nor; ++i) {
      ffMulRow(m1, FF_ZERO, noc);
      for (uint32_t j = 0; j < noc; ++j) {
         FEL a = readFel();
         ffInsert(m1, j, a);
      }
      ffWriteRows(out, m1, 1, noc);
      assertEndOfLine();
   }
   sysFree(m1);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

static void convertIntegerMatrix()
{
   uint32_t nor = MTX_NVAL;
   uint32_t noc = MTX_NVAL;
   while (1) {
      if (tryParseHeaderField(&nor, "rows=")) continue;
      if (tryParseHeaderField(&nor, "nor=")) continue;
      if (tryParseHeaderField(&noc, "cols=")) continue;
      if (tryParseHeaderField(&noc, "noc=")) continue;
      if (*lptr == 0) break;
      mtxAbort(MTX_HERE,"Unknown header field");
   }
   if (nor == MTX_NVAL) {mtxAbort(MTX_HERE, "Missing header field \"rows\".");}
   if (noc == MTX_NVAL) {mtxAbort(MTX_HERE, "Missing header field \"cols\".");}

   MESSAGE(0,("%dx%d integer matrix\n",nor,noc));
   int32_t* x = NALLOC(int32_t, noc);
   WriteHeader(MTX_TYPE_INTMATRIX, nor, noc);
   for (uint32_t i = 1; i <= nor; ++i) {
      for (uint32_t j = 0; j < noc; ++j) {
         x[j] = parse32s();
      }
      sysWrite32(out,x,noc);
      assertEndOfLine();
   }
   free(x);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

static void convertPermutation()
{
   uint32_t degree = MTX_NVAL;
   while (1) {
      if (tryParseHeaderField(&degree, "degree=")) continue;
      if (tryParseHeaderField(&degree, "deg=")) continue;
      if (*lptr == 0) break;
      mtxAbort(MTX_HERE,"Unknown header field");
   }
   if (degree == MTX_NVAL) {
      mtxAbort(MTX_HERE,"Missing header field \"degree\".");
   }

   MESSAGE(0,("Permutation on %lu points\n",(unsigned long) degree));
   uint32_t *buf = NALLOC(uint32_t,degree);
   WriteHeader(MTX_TYPE_PERMUTATION, degree, 1);
   for (uint32_t i = 0; i < degree; ++i) {
      uint32_t point = parse32u();
      if (point == 0 || point > degree) {
         mtxAbort(MTX_HERE,"Invalid point %lu in permutation of degree %lu",
               (unsigned long) point, (unsigned long) degree);
      }
      buf[i] = point - 1;
   }
   sysWrite32(out, buf, degree);
   sysFree(buf);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

static void convertPolynomial()
{
   uint32_t fieldOrder = MTX_NVAL;
   int32_t degree = -2;
   int hasDegree = 0;

   while (1) {
      if (tryParseHeaderField32s(&degree, "degree=")) {
         hasDegree = 1;
         continue;
      }
      if (tryParseHeaderField(&fieldOrder, "field=")) {continue;}
      if (*lptr == 0) {break;}
      mtxAbort(MTX_HERE, "Unknown header field");
   }
   if (!hasDegree) {mtxAbort(MTX_HERE, "Missing header field \"degree\".");}
   if (fieldOrder == MTX_NVAL) {mtxAbort(MTX_HERE, "Missing header field \"field\".");}

   MESSAGE(0, ("Polynomial of degree %ld over GF(%lu)\n",
               (long) degree, (unsigned long)fieldOrder));
   ffSetField(fieldOrder);
   Poly_t* p = polAlloc(fieldOrder, degree);
   for (int32_t i = 0; i <= degree; ++i) {
      p->data[i] = readFel();
   }
   polWrite(p, out);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

/// Convert one member.

static int convert(void)
{
   if (!tryReadLine())
      return 0;
   char *c;

   // Check for new file format
   if ((c = strstr(lbuf,"MeatAxeFileInfo")) != NULL) {
      char *d = lbuf;
      for (c += 15; *c != 0 && *c != '"'; ++c) {
      }
      if (*c != '"') {
         mtxAbort(MTX_HERE,"%s: Bad file format",inpname);
      }
      for (++c; *c != '"' && *c != 0; ++c) {
         *d++ = *c;
      }
      *d = 0;
      GrpLibFormat = 1;
   }

   if (tryParseLiteral("matrix")) {
      convertMatrix();
   } else if (tryParseLiteral("integer matrix") || tryParseLiteral("integer-matrix")) {
      convertIntegerMatrix();
   } else if (tryParseLiteral("permutation")) {
      convertPermutation();
   } else if (tryParseLiteral("polynomial")) {
      convertPolynomial();
   } else {
      mtxAbort(MTX_HERE,"%s: Unrecognized object header",inpname);
   }
   assertEndOfLine();
   return 1;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{
   App = appAlloc(&AppInfo,argc,argv);
   appGetArguments(App,2,2);

   // input file
   inpname = App->argV[0];
   if (strcmp(inpname,"-")) {
      src = sysFopen(inpname, "r::lib");
      if (src == NULL)
         mtxAbort(MTX_HERE,"Cannot open %s",inpname);
   } else {
      src = stdin;
   }

   // output file
   outname = App->argV[1];
   out = sysFopen(outname,"wb");
   if (out == NULL)
      mtxAbort(MTX_HERE,"Cannot open %s for output",outname);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   init(argc,argv);
   int inputFileScope = mtxBeginScope(provideInputFilePosition, NULL);
   while (convert()) {
      ++MemberCount;
   }
   if (MemberCount == 0) {
      MESSAGE(0,("Warning: %s is empty\n",inpname));
   }
   mtxEnd(inputFileScope);
   fclose(out);
   appFree(App);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// *INDENT-OFF*

/**
@page prog_zcv zcv - Convert Text to Binary Format

@see @ref prog_zpr

@section zcv_syntax Command Line
<pre>
zcv @em Options @em TextFile @em DataFile
</pre>

@par @em Options
Standard options, see @ref prog_stdopts

@par @em TextFile
Input file (text)

@par @em DataFile
Output file (binary)

@section zcv_inp Input Files

@par @em TextFile
Input file (text)

@section zcv_out Output Files

@par @em DataFile
Output file (binary)

@section zcv_desc Description
This program converts a text file into binary format.
If the input file name is "-", input is read from stdin.

@par Text File Format
The text file is interpreted line by line. Empty lines are ignored
completetly. If the line contains one or more "#", which may occur at
any position, the first "#" and all remaining characters on this line
are ingored.

Each object, for example a matrix or a permutation, consists of a one
line header followed by the data. Here are some examples of possible
header formats:
<pre>
matrix field=16 rows=10 cols=10
permutation degree=10026
polynomial field=2 degree=23
</pre>
The header may be given in a different format, for example
<pre>
MeatAxeFileInfo := "matrix field=5 rows=100 cols=100";
</pre>
*/

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
