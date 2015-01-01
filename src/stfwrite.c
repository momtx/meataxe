////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Structured Text File (STF) output functions.
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#define MaxCharsPerLine 80


MTX_DEFINE_FILE_INFO 

/// @addtogroup stf
/// @{


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a raw value.

int StfPut(StfData *f, const char *text)
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

int StfPutInt(StfData *f, int value)
{
    char tmp[20];
    sprintf(tmp,"%d",value);
    return StfPut(f,tmp);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a string.

int StfPutString(StfData *f, const char *text)
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
    result = StfPut(f,tmp);
    FREE(tmp);
    return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a vector.

int StfPutVector(StfData *f, int size, const int *value)
{
    int i;
    if (value == NULL || size < 0 || size > 100000)
	return -1;
    if (f == NULL || f->File == NULL)
	return -1;
    StfPut(f,"[");
    for (i = 0; i < size; ++i)
    {
	StfPutInt(f,value[i]);
	if (i < size - 1)
	    StfPut(f,",");
    }
    StfPut(f,"]");
    return 0;
}





////////////////////////////////////////////////////////////////////////////////////////////////////
/// Start a new entry.
/// This function begins a new entry. Be sure to terminate any incomplete 
/// entries with StfEndEntry() before you start a new entry. If you don't, 
/// the incomplete entry may be lost, and the data file may become corrupted.
///
/// Before you use %StfBeginEntry(), check if you can do the job with one of
/// the StfWriteXXX() functions. In particular, there are functions to write
/// integers and sequences of integers. If you have more complicated data to write,
/// you may need to construct the output manually.
/// Here is an example:
/// @code
/// StfBeginEntry(f,"Param");
/// StfPut(f,"(");
/// StfPutInt(f,11);
/// StfPut(f,":");
/// StfPutInt(f,22);
/// StfPut(f,")");
/// StfEndEntry(f);
/// @endcode
/// This code produces the following line in the output file:
/// <pre>
/// Param := (11:22);
/// </pre>
/// @param f Pointer to a structured text file (STF) object.
/// @param name Name of the entry.
/// @return 0 on success, -1 on error.

int StfBeginEntry(StfData *f, const char *name)
{
    if (name == NULL)
	return -1;
    if (f == NULL || f->File == NULL)
	return -1;
    if (StfPut(f,name) || StfPut(f," := "))
	return -1;
    return 0;
}




////////////////////////////////////////////////////////////////////////////////////////////////////
/// End entry.
/// This function terminates the current entry and flushes the STF's line
/// buffer. See StfBeginEntry() for an example.
/// @param f Pointer to a structured text file (STF) object.
/// @return 0 on success, -1 on error.

int StfEndEntry(StfData *f)
{
    if (f == NULL || f->File == NULL)
	return -1;
    return StfPut(f,";\n");
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a string.
/// This function writes an arbitrary text to a structured text file. 
/// For example, the statement
/// @code
/// StfWriteValue(f,"Note","This is a note");
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

int StfWriteValue(StfData *f, const char *name, const char *value)
{
    if (name == NULL || value == NULL)
	return -1;
    if (f == NULL || f->File == NULL)
	return -1;
    if (StfBeginEntry(f,name) != 0)
	return -1;
    StfPut(f,value);
    StfEndEntry(f);
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
/// StfWriteValue(f,"Title","This is a test ");
/// @endcode
/// produces the following output line:
/// <pre>
/// Title := "This is a test ";
/// </pre>
/// Unlike StrWriteValue() this function preserves leading and trailing spaces.

int StfWriteString(StfData *f, const char *name, const char *value)
{
    if (name == NULL || value == NULL)
    {
	MTX_ERROR("name or value invalid");
	return -1;
    }
    if (f == NULL || f->File == NULL)
    {
	MTX_ERROR("Invalid file");
	return -1;
    }
    if (StfBeginEntry(f,name) != 0)
	return -1;
    StfPutString(f,value);
    StfEndEntry(f);
    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write an enteger.
/// This function writes an integer to a structured text file. For example, 
/// the statement
/// @code
/// StfWriteInt(f,"Dimension",42);
/// @endcode
/// produces the following output line:
/// <pre>
/// Dimension := 42;
/// </pre>
/// @param f Pointer to a structured text file (STF) object.
/// @param name Name of the entry.
/// @param value Value of the entry.
/// @return 0 on success, -1 on error.

int StfWriteInt(StfData *f, const char *name, int value)
{
    if (name == NULL)
	return -1;
    if (f == NULL || f->File == NULL)
    {
	MTX_ERROR1("f: %E",MTX_ERR_BADARG);
	return -1;
    }
    if (StfBeginEntry(f,name) != 0)
	return -1;
    StfPutInt(f,value);
    StfEndEntry(f);
    return 0;
}




////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write a vector.
/// This function writes a sequence of integers to a structured text file.
/// For example, the statement
/// @code
/// int dims[5] = {11,22,33,44,55};
/// StfWriteVector(f,"Dimensions",dims,5);
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

int StfWriteVector(StfData *f, const char *name, int size, const int *value)
{
    if (name == NULL || value == NULL || size < 0 || size > 100000)
	return -1;
    if (f == NULL || f->File == NULL)
	return -1;
    if (StfBeginEntry(f,name) != 0)
	return -1;
    StfPutVector(f,size,value);
    StfEndEntry(f);
    return 0;
}

/// @}
