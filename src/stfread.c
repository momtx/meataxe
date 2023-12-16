////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Structured Text File (STF) input functions.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>


/// @addtogroup stf
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read next line.
/// This function reads a single text line into the STF object's internal 
/// line buffer and prepares the text for parsing with stfGetXXX() functions.
/// stfReadLine() strips comments and assembles multi-line texts into a 
/// single line. Thus, the application need not handle comments and multi-line
/// texts.
/// @param f Pointer to a structured text file (STF) object.
/// @return 0 on success, -1 on end-of-file or error.

int stfReadLine(StfData *f)
{
    char lbuf[6000];
    int ch;
    int len;
    int tlen = 0;

    f->getPtr = "";
    while (!feof(f->file))
    {
	lbuf[0] = 0;
	if (fgets(lbuf,sizeof(lbuf),f->file) == NULL)
	    break;
	++f->lineNo;
	len = strlen(lbuf);
	if (lbuf[len-1] == '\n') 
	    --len;
	if (len <= 0) continue;
	if (lbuf[0] == '#') continue;
	if (tlen + len >= f->lineBufSize)
	{
	    char *newbuf;
	    int newsize = tlen + len + 5;
	    newbuf = NREALLOC(f->lineBuf,char,newsize);
	    if (newbuf == NULL)
		return -1;
	    f->lineBuf = newbuf;
	    f->lineBufSize = newsize;
	}
	strcpy(f->lineBuf + tlen,lbuf);
	tlen += len;

	ch = fgetc(f->file);
	if (ch != '\t' && ch != EOF)
	{
	    ungetc(ch,f->file);
	    break;
	}
    }
    f->lineBuf[tlen] = 0;
    return f->lineBuf[0] != 0 ? 0 : -1;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Get entry name.
/// This function extracts the name part of internal line buffer and prepares 
/// the buffer for further parsing with stfGetXXX() functions.
/// On return,  <tt>f->GetPtr</tt> points to the first non-space character after
/// the ":=".
/// %stfGetName() can be called only after a line was successfully read with
/// stfReadLine(). It must be called before any of the stfGetXXX()
/// functions. See stfGetInt() for an example.
/// @param f Pointer to a structured text file (STF) object.
/// @return Name found in text line or |NULL| on error.

const char *stfGetName(StfData *f)
{
    char *c, *name;
    
    f->getPtr = NULL;

    /* Skip leading whitespace */
    for (c = f->lineBuf; *c != 0 && isspace(*c); ++c);
    if (*c == 0) return NULL;

    /* Parse name */
    name = c;
    while (*c != 0 && !isspace(*c)) ++c;
    if (*c == 0)
    {
	f->getPtr = c;
	return name;
    }
    else
	*c++ = 0;

    /* Skip " := " */
    *c++ = 0;
    while (*c != 0 && (isspace(*c) || *c == ':' || *c == '=')) ++c;
    f->getPtr = c;

    return name;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int stfGetULong(StfData *f, unsigned long *buf)
{
    char *c = f->getPtr;

    if (c == NULL)
	return -1;

    // Skip leading whitespace
    while (*c != 0 && isspace(*c)) ++c;

    // Parse number */
    if (!isdigit(*c)) 
	mtxAbort(MTX_HERE,"Invalid unsigned integer: %.20s",f->getPtr);
    *buf = 0;
    while (*c != 0 && isdigit(*c))
    {
	*buf = *buf * 10 + (unsigned char) (*c - '0');
	++c;
    }

    f->getPtr = c;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read an integer.
/// This function gets one integer from the current line and increments the
/// read pointer accordingly. Before this function is called, a line must 
/// have been read with stfReadLine() and prepared with stfGetName().
///
/// Reading starts at the current position. Any leading blanks are skipped.
/// On return, the new current position is the character following the last 
/// digit. If there is no integer to read, the current position is not 
/// changed, and the function returns -1.
///
/// Here is an example:
/// @code
/// StfFile *f = stfOpen("test","r");
/// int dim, degree, result = 0;
/// while (result == 0 && stfReadLine(f) == 0)
/// {
///     char *name = stfGetName(f);
///     if (!strcmp(name,"Dimension"))
///         result = stfGetInt(f,&dim);
///     else if (!strcmp(name,"Degree"))
///         result = stfGetInt(f,&degree);
/// }
/// @endcode
/// This code fragment opens a text file and reads two parameters, "Degree" and 
/// "Dimension", into the variables @c degree and @c dim, respectively.
/// @param f Pointer to a structured text file (STF) object.
/// @param buf Pointer to a buffer receiving the value.
/// @return 0 on success, -1 on error.

int stfGetInt(StfData *f, int *buf)
{
    char *c = f->getPtr;
    int neg = 0;

    if (c == NULL)
	return -1;

    /* Skip leading whitespace */
    while (*c != 0 && isspace(*c)) ++c;

    /* Parse sign */
    if (*c == '-') 
    { 
	neg = 1; 
	++c; 
    }

    /* Parse number */
    if (!isdigit(*c)) 
    {
	mtxAbort(MTX_HERE,"Invalid integer: %.20s",f->getPtr);
	return -1;
    }
    *buf = 0;
    while (*c != 0 && isdigit(*c))
    {
	*buf = *buf * 10 + (unsigned char) (*c - '0');
	++c;
    }
    if (neg)
	*buf = -*buf;

    /* Done */
    f->getPtr = c;
    return 0;
}


int stfGetU32(StfData *f, uint32_t *buf)
{
   unsigned long value;
   if (stfGetULong(f, &value) != 0)
      return -1;
   return (uint32_t) value;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read a string.
/// This function gets a string from the current line and increments the
/// read pointer accordingly. Before this function is called, a line must 
/// have been read with stfReadLine() and prepared with stfGetName().
/// The string is expected at the current position of the test file and must
/// be in C syntax, i.e., enclosed in double quotation marks.
/// @param f Pointer to a structured text file (STF) object.
/// @param buf Pointer to a buffer receiving the value.
/// @param bufsize Buffer size in bytes.
/// @return 0 on success, -1 on error.

int stfGetString(StfData *f, char *buf, size_t bufsize)
{
    char *c, *d, *s;

    /* Find start of string.
       --------------------- */
    for (c = f->getPtr; *c != 0 && *c != '"' && isspace(*c); ++c)
	;
    if (*c != '"')
    {
	mtxAbort(MTX_HERE,"Missing \"");
	return -1;
    }
    s = c;

    /* Traverse the string, replacing escape sequences.
       ------------------------------------------------ */
    for (d = c + 1; *d != 0 && *d != '"'; )
    {
	if (*d == '\\')
	{
	    ++d;
	    switch (*d)
	    {
		case 'n': *c++ = '\n'; break;
		case 'r': *c++ = '\r'; break;
		case 't': *c++ = '\t'; break;
		case 'a': *c++ = '\a'; break;
		case 'b': *c++ = '\b'; break;
		case 'f': *c++ = '\f'; break;
		case '"': *c++ = '"'; break;
		default:
		    mtxAbort(MTX_HERE,"Line %d: Invalid escape sequence in string",
			f->lineNo);
		    return -1;
	    }
	    ++d;
	}
	else 
	    *c++ = *d++;
    }
    if (*d != '"')
    {
	mtxAbort(MTX_HERE,"Line %d: Unexpected end of line in string",f->lineNo);
	return -1;
    }

    /* Check buffer size.
       ------------------ */
    if (bufsize < c - s + 1)
    {
	mtxAbort(MTX_HERE,"Line %d: Buffer too small",f->lineNo);
	return -1;
    }

    /* Copy the string.
       ---------------- */
    memcpy(buf,s,c - s);
    buf[c - s] = 0;

    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Skip text.
/// This function reads (and skips) the text given by @a pattern.
/// Before using this function, a line must have been read with stfReadLine()
/// and prepared with stfGetname(). Reading starts at the current position, 
/// i.e., at the first non-space character after the ":=", or at the first
/// character that was not read by the previous stfGetXXX() or stfMatch().
/// A space in @a pattern matches any number (including 0) of spaces or tabs. 
/// Any other characters in @a pattern are matched one-to-one against the input 
/// line. 
/// If @a pattern is matched completely, the current position is updated 
/// to the character after the last matched character. Otherwise, the current 
/// position is not changed and the function returns -1.
/// @param f Pointer to a structured text file (STF) object.
/// @param pattern The text to be skipped.
/// @return 0, if the complete text in @a pattern has beed
/// skipped. -1 otherwise.

int stfMatch(StfData *f, const char *pattern)
{
    char *b = f->getPtr;

    if (b == NULL)
	return -1;

    for ( ; *b != 0 && *pattern != 0; ++pattern)
    {
	if (*pattern == ' ')
	{
	    while (*b != 0 && isspace(*b)) 
		++b;
	}
	else
	{
	    if (*b == *pattern) 
		++b;
	    else 
		return -1;
	}
    }
    if (*pattern == 0)
    {	
	f->getPtr = b;
	return 0;
    }
    return -1;
}




////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read a vector.
/// This function reads a sequence of integers. The sequence must have been
/// written with stfWriteVector() or at least in the same format.
///
/// Before using this function, a line must have been read with stfReadLine()
/// and prepared with stfGetname(). Reading starts at the current position, 
/// i.e., at the first non-space character after the ":=", or at the first
/// character that was not read by the previous stfGetXXX() or stfMatch().
///
/// The caller must supply two buffers, the data buffer (@a buf) and an 
/// integer buffer (@a bufsize). When %stfGetVector() is called, the variable
/// pointed to by @a bufsize must contain the maximal number of integers
/// that can be stored in @a buf. On successful return, the variable contains 
/// the number of integers that were actually stored, which may be smaller than
/// the original value.
/// If the vector is to long to fit into the user-supplied buffer, the
/// function reads as many entries as possible and returns -1.
///
/// Here is an example:
/// @code
/// char *name = stfGetName(f);
/// if (!strcmp(name,"Vector"))
/// {
///     int vec[10];
///     int vecsize = 10;
///     stfGetVector(f,&vecsize,vec);
///     printf("%d elements read\n",*vecsize);
/// }
/// @endcode
/// @param f Pointer to a structured text file (STF) object.
/// @param bufsize Pointer to a variable containing the buffer size. On return, the variable
///      contains the number of elements actually read.
/// @param buf Pointer to a buffer receiving the data.
/// @return The function returns $0$ on success and $-1$ on error.

int stfGetVector(StfData *f, int *bufsize, int *buf)
{
    char *c = f->getPtr;
    int i;
    if (buf == NULL || bufsize == NULL || *bufsize <= 0)
	return -1;
    
    if (stfMatch(f," ["))
	return (f->getPtr = c, -1);
    for (i = 0; i < *bufsize; ++i, ++buf)
    {
	if (stfMatch(f," ]") == 0)
	    break;
	if (i > 0 && stfMatch(f,","))
	    return (f->getPtr = c, -1);
	if (stfGetInt(f,buf))
	    break;
    }
    if (i >= *bufsize && stfMatch(f,"]"))
	return (f->getPtr = c, -1);
    *bufsize = i;
    return 0;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
