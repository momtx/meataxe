////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Convert a matrix or permutation from ASCII (readable)
////////////////////////////////////////////////////////////////////////////////////////////////////

#define MAXLINE 4000    // max. input line size

#include "meataxe.h"
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Global data


static int GrpLibFormat = 0;    /* File is in Group Library Format */
static uint32_t fl;
static int mod;
static uint32_t nor;
static uint32_t noc;
static FILE *src = NULL;                /* Input file */
static FILE *out;                       /* Output file */
static char lbuf[MAXLINE] = {0};        /* Input line buffer */
static char *lptr = lbuf;               /* Read pointer */
static const char *inpname = "[stdin]";
static const char *outname = "";
static int MemberCount = 0;             /* Number of members */

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

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read next input line, skip empty lines and strip comments.
/// Returns 0 on success, 1 on end-of-file.

static int readline()
{
   char *c;
   int mt = 1;

   while (mt) {
      lbuf[0] = 0;
      if (feof(src)) { return 1; }
      fgets(lbuf,sizeof(lbuf),src);
      if (ferror(src)) {
         mtxAbort(MTX_HERE,"Unexpected end of input file");
         return -1;
      }
      for (c = lbuf; *c != 0 && *c != '#'; ++c) {
         if (!isspace(*c)) { mt = 0; }}
      *c = 0;
   }
   lptr = lbuf;
   return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read integer number

static int32_t readlong()
{
   long l;
   int minus = 0;

   while (!isdigit(*lptr) && *lptr != '-') {
      while (*lptr != 0 && !isdigit(*lptr) && *lptr != '-') {
         ++lptr;
      }
      if (*lptr == 0) {
         if (readline()) {
            mtxAbort(MTX_HERE,"Unexpected end of input file");
            return -1;
         }
      }
   }
   if (*lptr == '-') { minus = 1; ++lptr; }
   if (!isdigit(*lptr)) {
      mtxAbort(MTX_HERE,"%s: Bad file format",inpname);
      return -1;
   }
   for (l = 0; isdigit(*lptr); ) {
      l *= 10;
      switch (*lptr) {
         case '0': break;
         case '1': l += 1;
            break;
         case '2': l += 2;
            break;
         case '3': l += 3;
            break;
         case '4': l += 4;
            break;
         case '5': l += 5;
            break;
         case '6': l += 6;
            break;
         case '7': l += 7;
            break;
         case '8': l += 8;
            break;
         case '9': l += 9;
            break;
      }
      ++lptr;
      if ((mod != 2) && (mod != 5) && (fl >= 2) && (fl <= 9)) { break; }
   }
   return minus ? -l : l;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void WriteHeader(uint32_t a, uint32_t b, uint32_t c)
{
   uint32_t hdr[3] = {a, b, c};
   sysWrite32(out,hdr,3);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static const char* strPrefix(const char* s, const char* prefix)
{
   if (s == NULL || prefix == NULL)
      return NULL;
   while (*prefix != 0 && *s == *prefix) {
      ++s;
      ++prefix;
   }
   return (*prefix == 0) ? s : NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static int parse32u(uint32_t* var, const char *s, const char* prefix)
{
   const char* valueStr = strPrefix(s, prefix);
   if (valueStr == NULL) return 0;
      
   char* end = NULL;
   unsigned long value = strtoul(valueStr, &end, 10);
   if (*valueStr != 0 && *end == 0) {
      if (value <= 0xFFFFFFFFU) {
         *var = (uint32_t) value;
         return 1;
      }
   }
   
   mtxAbort(MTX_HERE, "%s: invalid number in object header: \"%s\"", inpname, s);
   return 0;
}
         
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Convert matrix

static void ConvertMatrix()
{
   long i, j;
   PTR m1;

   MESSAGE(0,("%dx%d matrix over GF(%d)\n",nor,noc,fl));
   ffSetField(fl);
   m1 = ffAlloc(1, noc);
   WriteHeader(fl,nor,noc);
   for (i = 1; i <= nor; ++i) {
      ffMulRow(m1,FF_ZERO, noc);
      for (j = 0; j < noc; ++j) {
         long val = readlong();
         ffInsert(m1,j,ffFromInt(val));
      }
      ffWriteRows(out, m1, 1, noc);
   }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Convert integer matrix

static void ConvertIntMatrix()
{
   long i, j;
   int32_t *x;

   MESSAGE(0,("%dx%d integer matrix\n",nor,noc));
   x = NALLOC(int32_t,noc);
   WriteHeader(MTX_TYPE_INTMATRIX, nor, noc);         /* 8 = T_IMAT */
   for (i = 1; i <= nor; ++i) {
      for (j = 0; j < noc; ++j) {
         x[j] = readlong();
      }
      sysWrite32(out,x,noc);
   }
   free(x);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

/// Convert permutation (new format)

static void ConvertPermutation()
{
   long i;
   uint32_t *buf;
   long kk;

   MESSAGE(0,("Permutation on %d points\n",nor));
   buf = NALLOC(uint32_t,nor);
   WriteHeader(MTX_TYPE_PERMUTATION, nor, 1);

   for (i = 0; i < nor; ++i) {
      kk = readlong();
      MTX_ASSERT(kk > 0);
      MTX_ASSERT(kk <= nor);
      buf[i] = kk - 1;
   }
   sysWrite32(out, buf, nor);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Convert polynomial (new format)

static void ConvertPolynomial()
{
   long i;
   Poly_t *p;

   MESSAGE(0,("Polynomial of degree %d over GF(%d)\n",nor,fl));
   p = polAlloc(fl,nor);
   out = sysFopen(outname,"wb");
   for (i = 0; i <= nor; ++i) {
      long kk = readlong();
      p->Data[i] = ffFromInt(kk);
   }
   polWrite(p,out);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
///  Convert permutation (mode 2)

void convperm()
{
   long i, val;
   PTR m1;

   MESSAGE(0,("%dx%d permutation matrix over GF(%d)\n",nor,noc,fl));
   ffSetField(fl);
   m1 = ffAlloc(1, noc);
   WriteHeader(fl,nor,noc);
   for (i = 1; i <= nor; ++i) {
      val = readlong();
      ffMulRow(m1,FF_ZERO, noc);
      ffInsert(m1,val - 1,FF_ONE);
      ffWriteRows(out, m1, 1, noc);
   }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read the next header. Returns 0 on end of file, 1 on success.

static int ReadHeader(void)
{
   if (readline()) {
      return 0;
   }
   return 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

/// Convert one member.

static void Convert(void)
{
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

   if (!strncmp(lbuf,"matrix",6)) {
      fl = nor = noc = MTX_NVAL;
      for (char* c = strtok(lbuf + 6," \t\n"); c != NULL; c = strtok(NULL," \t\n")) {
         if (parse32u(&fl, c, "field=")) continue;
         if (parse32u(&nor, c, "nor=")) continue;
         if (parse32u(&nor, c, "rows=")) continue;
         if (parse32u(&noc, c, "noc=")) continue;
         if (parse32u(&noc, c, "cols=")) continue;
         mtxAbort(MTX_HERE,"%s: Bad object header: \"%s\"", inpname, c);
      }
      if ((nor == MTX_NVAL)
            || (noc == MTX_NVAL)
            || (fl == MTX_NVAL)) {
         mtxAbort(MTX_HERE,"%s: Bad file format (missing header field)",inpname);
      }
      readline();
      ConvertMatrix();
      return;
   }
   if (!strncmp(lbuf,"integer matrix",14) || !strncmp(lbuf,"integer-matrix",14)) {
      char *c;
      nor = noc = MTX_NVAL;
      for (c = strtok(lbuf + 14," \t\n"); c != NULL; c = strtok(NULL," \t\n")) {
         if (parse32u(&nor, c, "nor=")) continue;
         if (parse32u(&nor, c, "rows=")) continue;
         if (parse32u(&noc, c, "noc=")) continue;
         if (parse32u(&noc, c, "cols=")) continue;
         mtxAbort(MTX_HERE,"%s: Bad object header: \"%s\"", inpname, c);
      }
      if (nor == MTX_NVAL || noc == MTX_NVAL) {
         mtxAbort(MTX_HERE,"%s: Bad file format (missing header field)",inpname);
      }
      readline();
      ConvertIntMatrix();
      return;
   } else if (!strncmp(lbuf,"permutation",11)) {
      char *c;
      fl = nor = MTX_NVAL;
      noc = 1;
      for (c = strtok(lbuf + 11," \t\n"); c != NULL; c = strtok(NULL," \t\n")) {
         if (parse32u(&nor, c, "degree=")) continue;
         if (parse32u(&nor, c, "deg=")) continue;
         mtxAbort(MTX_HERE,"%s: Bad object header: \"%s\"", inpname, c);
      }
      if (nor == MTX_NVAL) {
         mtxAbort(MTX_HERE,"%s: Bad file format (missing header field)",inpname);
      }
      readline();
      ConvertPermutation();
      return;
   } else if (!strncmp(lbuf,"polynomial",10)) {
      char *c;
      fl = nor = MTX_NVAL;
      noc = 1;
      for (c = strtok(lbuf + 11," \t\n"); c != NULL; c = strtok(NULL," \t\n")) {
         if (parse32u(&nor, c, "degree=")) continue;
         if (parse32u(&nor, c, "deg=")) continue;
         if (parse32u(&fl, c, "field=")) continue;
         mtxAbort(MTX_HERE,"%s: Bad object header: \"%s\"", inpname, c);
      }
      if (fl == MTX_NVAL || noc == MTX_NVAL) {
         mtxAbort(MTX_HERE,"%s: Bad file format (missing header field)",inpname);
      }
      readline();
      ConvertPolynomial();
      return;
   } else {
      mtxAbort(MTX_HERE,"%s: Unsupported file format",inpname);
   }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////

static int Init(int argc, char **argv)
{
   if ((App = appAlloc(&AppInfo,argc,argv)) == NULL) {
      return -1;
   }
   if (appGetArguments(App,2,2) < 0) {
      return -1;
   }

   // input file
   inpname = App->ArgV[0];
   if (strcmp(inpname,"-")) {
      src = sysFopen(inpname, "r::lib");
      if (src == NULL) {
         mtxAbort(MTX_HERE,"Cannot open %s",inpname);
         return -1;
      }
   } else {
      src = stdin;
   }

   // output file
   outname = App->ArgV[1];
   out = sysFopen(outname,"wb");
   if (out == NULL) {
      mtxAbort(MTX_HERE,"Cannot open %s for output",outname);
      return -1;
   }

   return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////

static void Cleanup()
{
   fclose(out);
   if (App != NULL) { appFree(App); }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   if (Init(argc,argv) != 0) {
      return 0;
   }

   int context = mtxBegin("Converting %s", inpname);
   while (ReadHeader()) {
      Convert();
      ++MemberCount;
   }
   if (MemberCount == 0) {
      MESSAGE(0,("Warning: %s is empty",inpname));
   }
   mtxEnd(context);

   Cleanup();
   return 0;
}


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
