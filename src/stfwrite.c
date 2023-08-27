////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Structured Text File (STF) output functions.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#define MaxCharsPerLine 80


/// @addtogroup stf
/// @{


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a raw value.

int stfPut(StfData *f, const char *text)
{
    int len = strlen(text);
    if (len == 0)
	return 0;

    if ((f->OutPos + len) > MaxCharsPerLine)
    {
	fputs("\n\t",f->File);
	f->OutPos = 8;
	++f->LineNo;
    }
    fputs(text,f->File);
    f->OutPos += len;
    if (text[len-1] == '\n')
    {
	f->OutPos = 0;
	++f->LineNo;
    }
    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write an integer.

int stfPutInt(StfData *f, int value)
{
    char tmp[20];
    sprintf(tmp,"%d",value);
    return stfPut(f,tmp);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a string.

int stfPutString(StfData *f, const char *text)
{
    char *tmp = NALLOC(char,2 * strlen(text) + 5);
    char *c = tmp;
    const char *t = text;
    int result;

    *c++ = '"';
    while (*t != 0)
    {
	switch (*t)
	{
	    case '\n': *c++ = '\\'; *c++ = 'n'; break;
	    case '\r': *c++ = '\\'; *c++ = 'r'; break;
	    case '\t': *c++ = '\\'; *c++ = 't'; break;
	    case '\a': *c++ = '\\'; *c++ = 'a'; break;
	    case '\b': *c++ = '\\'; *c++ = 'b'; break;
	    case '\f': *c++ = '\\'; *c++ = 'f'; break;
	    case '"': *c++ = '\\'; *c++ = '"'; break;
	    default:
		*c++ = *t;
	}
	++t;
    }
    *c++ = '"';
    *c = 0;
    result = stfPut(f,tmp);
    sysFree(tmp);
    return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a vector.

int stfPutVector(StfData *f, int size, const int *value)
{
    int i;
    if (value == NULL || size < 0 || size > 100000)
	return -1;
    if (f == NULL || f->File == NULL)
	return -1;
    stfPut(f,"[");
    for (i = 0; i < size; ++i)
    {
	stfPutInt(f,value[i]);
	if (i < size - 1)
	    stfPut(f,",");
    }
    stfPut(f,"]");
    return 0;
}





////////////////////////////////////////////////////////////////////////////////////////////////////
/// Start a new entry.
/// This function begins a new entry. Be sure to terminate any incomplete 
/// entries with stfEndEntry() before you start a new entry. If you don't, 
/// the incomplete entry may be lost, and the data file may become corrupted.
///
/// Before you use %stfBeginEntry(), check if you can do the job with one of
/// the stfWriteXXX() functions. In particular, there are functions to write
/// integers and sequences of integers. If you have more complicated data to write,
/// you may need to construct the output manually.
/// Here is an example:
/// @code
/// stfBeginEntry(f,"Param");
/// stfPut(f,"(");
/// stfPutInt(f,11);
/// stfPut(f,":");
/// stfPutInt(f,22);
/// stfPut(f,")");
/// stfEndEntry(f);
/// @endcode
/// This code produces the following line in the output file:
/// <pre>
/// Param := (11:22);
/// </pre>
/// @param f Pointer to a structured text file (STF) object.
/// @param name Name of the entry.
/// @return 0 on success, -1 on error.

int stfBeginEntry(StfData *f, const char *name)
{
    if (name == NULL)
	return -1;
    if (f == NULL || f->File == NULL)
	return -1;
    if (stfPut(f,name) || stfPut(f," := "))
	return -1;
    return 0;
}




////////////////////////////////////////////////////////////////////////////////////////////////////
/// End entry.
/// This function terminates the current entry and flushes the STF's line
/// buffer. See stfBeginEntry() for an example.
/// @param f Pointer to a structured text file (STF) object.
/// @return 0 on success, -1 on error.

int stfEndEntry(StfData *f)
{
    if (f == NULL || f->File == NULL)
	return -1;
    return stfPut(f,";\n");
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a string.
/// This function writes an arbitrary text to a structured text file. 
/// For example, the statement
/// @code
/// stfWriteValue(f,"Note","This is a note");
/// @endcode
/// produces the following output line:
/// <pre>
/// Note := This is a note;
/// </pre>
/// Note that any leading spaces in the value will be stripped off when reading the file.
/// @param f Pointer to a structured text file (STF) object.
/// @param name Name of the entry.
/// @param value Value of the entry.
/// @return 0 on success, -1 on error.

int stfWriteValue(StfData *f, const char *name, const char *value)
{
    if (name == NULL || value == NULL)
	return -1;
    if (f == NULL || f->File == NULL)
	return -1;
    if (stfBeginEntry(f,name) != 0)
	return -1;
    stfPut(f,value);
    stfEndEntry(f);
    return 0;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a string.
/// @param f Pointer to a structured text file (STF) object.
/// @param name Name of the entry.
/// @param value String to write.
/// @return 0 on success, -1 on error.
/// This function writes an arbitrary string to a structured text file. 
/// For example, the statement
/// @code
/// stfWriteValue(f,"Title","This is a test ");
/// @endcode
/// produces the following output line:
/// <pre>
/// Title := "This is a test ";
/// </pre>
/// Unlike StrWriteValue() this function preserves leading and trailing spaces.

int stfWriteString(StfData *f, const char *name, const char *value)
{
    if (name == NULL || value == NULL)
    {
	mtxAbort(MTX_HERE,"name or value invalid");
	return -1;
    }
    if (f == NULL || f->File == NULL)
    {
	mtxAbort(MTX_HERE,"Invalid file");
	return -1;
    }
    if (stfBeginEntry(f,name) != 0)
	return -1;
    stfPutString(f,value);
    stfEndEntry(f);
    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write an enteger.
/// This function writes an integer to a structured text file. For example, 
/// the statement
/// @code
/// stfWriteInt(f,"Dimension",42);
/// @endcode
/// produces the following output line:
/// <pre>
/// Dimension := 42;
/// </pre>
/// @param f Pointer to a structured text file (STF) object.
/// @param name Name of the entry.
/// @param value Value of the entry.
/// @return 0 on success, -1 on error.

int stfWriteInt(StfData *f, const char *name, int value)
{
    if (name == NULL)
	return -1;
    if (f == NULL || f->File == NULL)
    {
	mtxAbort(MTX_HERE,"f: %s",MTX_ERR_BADARG);
	return -1;
    }
    if (stfBeginEntry(f,name) != 0)
	return -1;
    stfPutInt(f,value);
    stfEndEntry(f);
    return 0;
}




////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a vector.
/// This function writes a sequence of integers to a structured text file.
/// For example, the statement
/// @code
/// int dims[5] = {11,22,33,44,55};
/// stfWriteVector(f,"Dimensions",dims,5);
/// @endcode
/// produces the following output line:
/// <pre>
/// Dimensions := [11,22,33,44,55];
/// </pre>
/// @param f Pointer to a structured text file (STF) object.
/// @param name Name of the entry.
/// @param size Size of the vector (number of entries, not bytes).
/// @param value Pointer to the vector.
/// @return The function returns $0$ on success and $-1$ on error.

int stfWriteVector(StfData *f, const char *name, int size, const int *value)
{
    if (name == NULL || value == NULL || size < 0 || size > 100000)
	return -1;
    if (f == NULL || f->File == NULL)
	return -1;
    if (stfBeginEntry(f,name) != 0)
	return -1;
    stfPutVector(f,size,value);
    stfEndEntry(f);
    return 0;
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
