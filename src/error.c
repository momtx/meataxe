/* ============================= C MeatAxe ==================================
   File:        $Id: error.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Error handling
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>



static void (*ErrorHandler)(const MtxErrorRecord_t *err) = NULL;
static FILE *LogFile = NULL;



static void DefaultHandler(const MtxErrorRecord_t *e)
{
    static int count = 0;
    if (LogFile == NULL)
	LogFile = stderr;

    if (e->FileInfo != NULL)
	fprintf(LogFile,"%s(%d):",e->FileInfo->BaseName,e->LineNo);
    fprintf(LogFile,"%s\n",e->Text);
    if (--count <= 0)
	exit(255);
}


/**
!section general.err
 ** Define an application error handler.
!synopsis 
    typedef struct { ... } MtxErrorRecord;
    typedef void MtxErrorHandler_t(const MtxErrorRecord *err);
    MtxErrorHandler_t *MtxSetErrorHandler(MtxErrorHandler_t *h);
 ** @param h
    Address of the new error handler or |NULL| to restore the default handler.
 ** @return
    Returns the previous error handler.
!description
    This function defines an application error handler which is called every
    time an error occurs inside the {\MeatAxe} library.
 ** @see @errmsg MtxErrorRecord_t MtxErrorHandler_t
 **/

MtxErrorHandler_t *MtxSetErrorHandler(MtxErrorHandler_t *h)

{
    MtxErrorHandler_t *old = ErrorHandler;
    ErrorHandler = h;
    return old;
}



/**
 ** Signal an error.
 ** @param fi
    Pointer to a file information structure.
 ** @param line
    Line number where the error occured.
 ** @param text
    Error description.
 ** @param ...
    Optional arguments to be inserted into the error description.
 ** @return
    The return value is always 0.
!description
    This function is mainly used internally by the {\MeatAxe} library to
    report errors. If an application error handler has been installed, it 
    is executed. Otherwise an error message is written to the error log 
    file (stderr by default) and the program is terminated.

    The error message, |msg|, is formatted with |MtxFormatMessage()|.
 ** @see @errmsg
 **/

int MtxError(MtxFileInfo_t *fi, int line, const char *text, ...)

{
    va_list al;
    MtxErrorRecord_t err;
    char txtbuf[2000];

    /* Remove directory names from file name
       ------------------------------------- */
    if (fi != NULL && fi->BaseName == NULL)
    {
	const char *fn;
	for (fn = fi->Name; *fn != 0; ++fn);
	while (fn != fi->Name && fn[-1] != '/' && fn[-1] != '\\')
	    --fn;
	fi->BaseName = fn;
    }

    /* Create error record
       ------------------- */
    err.FileInfo = fi;
    err.LineNo = line;
    err.Text = txtbuf;
    va_start(al,text);
    MtxFormatMessage(txtbuf,sizeof(txtbuf),text,al);
    va_end(al);

    /* Call the error handler
       ---------------------- */
    if (ErrorHandler != NULL)
	ErrorHandler(&err);
    else
	DefaultHandler(&err);

    return 0;
}







/**
!structure MtxFileInfo_t "File information structure"
!synopsis 
typedef struct 
{ 
    const char *Name; 
    const char *BaseName; 
} MtxFileInfo_t;

!description
This data structure contains information about a source file. It is information
used when an error is reported via |MtxError()|. Each source file which uses the
{\MeatAxe} error mechanism, i.e.,|MtxError()| or any of the |MTX_ERROR| macros,
must also define one |MtxFileInfo_t| structure at file scope. There is a macro,
|MTX_DEFINE_FILE_INFO| which can be used to define this structure.
 ** @see MtxError MtxErrorRecord_t MTX_DEFINE_FILE_INFO
 **/


/**
!structure MtxErrorRecord_t "Error data structure"
!synopsis 
typedef struct 
{ 
    const MtxFileInfo_t *fi; 
    int LineNo; 
    const char *Text; 
} MtxErrorRecord_t;

!description
    This data structure contains detailed information on an error that occured
    inside the {\MeatAxe} library.
 ** @see MtxSetErrorHandler
 **/

/**
!structure MtxErrorHandler_t "Error handler function type"
!synopsis 
    typedef void MtxErrorHandler_t(const MtxErrorRecord_t *err);
!description
    This is the type of an application error handler.
 ** @see MtxSetErrorHandler MtxErrorRecord_t
 **/

/**
 ** Define file information.
!synopsis 
    MTX_DEFINE_FILE_INFO
!description
    This macro must be included in each source file that uses the
    {\MeatAxe} error reporting mechanism. It defines a data structure
    of type |MtxFileInfo_t| that is used by the error reporting macros
    |MTX_ERRORx|.
 ** @see MtxError
 **/

/**
!definesection general.err
Errors can occur inside the {\MeatAxe} library for many reasons:
\begin{itemize}
\item insufficient system resources (memory, disk space),
\item device or operating system failure,
\item invalid parameters passed to a {\MeatAxe} library function,
\item invalid usage of {\MeatAxe} library functions, for example missing
    initialization, etc.
\end{itemize}
Note that invalid parameters are not always detected by the {\MeatAxe} library.
For example, most kernel functions such as |FfAdd()| do not check their 
arguments for the sake of performance. Thus, calling these function with 
invalid arguments may produce random results or even crash the program. 
On the other hand, most higher-level functions like |MatAdd()| do some
plausibility checks on their arguments before processing them.

When an error is detected, the default action is to terminate the program
immediately with an appropriate error message. While this minimizes the chance 
of not noticing an error, it may be undesirable in some situations. For example,
the application may be able to recover from the error. In order to prevent 
the {\MeatAxe} library from terminating the program, the application can
define an application error handler. This is a function which is called on
each error.

To use the {\MeatAxe} error reporting mechanism, the following steps are
necessary:
\begin{enumerate}
\item Insert the macro |MTX_DEFINE_FILE_INFO| somewhere at the top of each 
    source file. This macro defines a |MtxFileInfo_t| structure which is used
    in reporting errors.
\item To report an error, use any of the |MTX_ERRORx| macros. There are several
    variants for different numbers of arguments. Note that the error message
    passed to |MTX_ERRORx()| is in printf style, but only a resticted set of 
    format specifiers is available. See |MtxFormatMessage()| for details.
\end{enumerate}
Here is a short example:
\begin{verbatim}
#include "meataxe.h"

MTX_DEFINE_FILE_INFO

int divide(int a, int b)
{
    if (b == 0)
    {
	MTX_ERROR("Division by 0");
	return 0;
    }
    return a / b;
}
\end{verbatim}
Note that you must not assume that |MTX_ERROR| terminates the program. In fact
|MTX_ERROR| may return if a user-defined error handler has been installed.
 **/
