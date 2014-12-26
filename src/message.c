/* ============================= C MeatAxe ==================================
   File:        $Id: message.c,v 1.2 2007-09-09 21:38:11 mringe Exp $
   Comment:     Messages.
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

/// @addtogroup app
/// @{


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

/// Message level.
/// This value determines which messages generated  by MSG0()...MSG4() are actually
/// output and which are suppressed.
int MtxMessageLevel = 0;

/// Last error code.
int MtxErrNo = 0;


/// Internal data structure.
static struct msg_struct { int ErrNo; char *smsg; } msgs[] =
{
    { MTX_ERR_NOMEM, "Not enough memory" },
    { MTX_ERR_GAME_OVER, "Time limit exceeded" },
    { MTX_ERR_DIV0, "Division by zero" },
    { MTX_ERR_FILEFMT, "Bad file format" },
    { MTX_ERR_BADARG, "Bad argument" },
    { MTX_ERR_RANGE, "Argument out of range" },

    { MTX_ERR_NOTECH, "Matrix not in echelon form" },
    { MTX_ERR_NOTSQUARE, "Matrix not square" },
    { MTX_ERR_INCOMPAT, "Incompatible objects" },

    { MTX_ERR_BADUSAGE, "Bad syntax, try `-help'" },
    { MTX_ERR_OPTION, "Bad usage of option, try `-help'" },
    { MTX_ERR_NARGS, "Bad number of arguments, try `-help'" },

    { MTX_ERR_NOTMATRIX, "Not a matrix" },
    { MTX_ERR_NOTPERM, "Not a permutation" },
    
    { -1, NULL }
};


/* ------------------------------------------------------------------
   errmsg()
   ------------------------------------------------------------------ */

static char *errmsg(int en)
{
    static char buf[30];
    struct msg_struct *x;

    for (x = msgs; x->ErrNo >= 0; ++x)
        if (x->ErrNo == en)
	    return x->smsg;

    sprintf(buf,"Unknown error %d",en);
    return buf;
}




/* ------------------------------------------------------------------
   FormatString() - Copy string to message buffer

   Description:
       This function is used internally by MtxFormatMessage().
   ------------------------------------------------------------------ */

static void FormatString(char **buf, int *count, int bufsize, const char *c)
{
    while (*c != 0 && *count < bufsize)
    {
	*(*buf)++ = *c++;
	++*count;
    }
}


/* ------------------------------------------------------------------
   FormatInteger() - Write decimal integer to message buffer

   Description:
       This function is used internally by MtxFormatMessage().
   ------------------------------------------------------------------ */

static void FormatInteger(char **buf, int *count, int bufsize, int val)
{
    char tmp[20];
    sprintf(tmp,"%d",val);
    FormatString(buf,count,bufsize,tmp);
}

static void FormatLong(char **buf, int *count, int bufsize, long int val)
{
    char tmp[20];
    sprintf(tmp,"%ld",val);
    FormatString(buf,count,bufsize,tmp);
}




/// Format a message.
/// This function formats a message using a |printf|-like syntax.
/// Only the follwing four format codes may be used within |msg|:
/// - @c @%d prints a signed decimal integer. The corresponding argument
///   must be of type @c int.
/// - @c @%s prints a string. The corresponding argument must be a pointer
///   to a null-terminated string.
/// - @c @%E takes an @c int argument, which must be one the MeatAxe error
///   codes (@c MTX_ERR_xxxx) defined in "meataxe.h". It prints a 
///   description of the error.
/// - @c @%S prints the system error name corresponding to the current
///   value of @c errno.
///
/// @param buf Buffer for the message.
/// @param bufsize Size of buffer.
/// @param msg The message text.
/// @param al Optional arguments for the message (see description).

int MtxFormatMessage(char *buf, int bufsize, const char *msg, va_list al)
{
    char *d;
    int count = 0;

    d = buf;
    while (*msg != 0 && count < bufsize - 1)
    {
	if (*msg == '%')
	{
	    switch (*++msg)
	    {
		case 'l':
		    FormatLong(&d,&count,bufsize,va_arg(al,long));
		    break;
		case 'd':
		    FormatInteger(&d,&count,bufsize,va_arg(al,int));
		    break;
		case 's':
		    FormatString(&d,&count,bufsize,va_arg(al,char *));
		    break;
		case 'E':
		    FormatString(&d,&count,bufsize,errmsg(va_arg(al,int)));
		    break;
		case 'S':
		    FormatString(&d,&count,bufsize,strerror(errno));
		    break;
	    }
	    ++msg;
	}
	else
	{
	    if (count < bufsize - 1)
		*d++ = *msg++;
	}
    }

    /* Terminate message and return number of characters written
       --------------------------------------------------------- */
    *d = 0;
    return count;
}



/// Print a message.
/// This function writes a message to a file using a @c printf()-like syntax.
/// See MtxFormatMessage() for details.
/// @param f File to write to.
/// @param fmt Message to write
/// @param ...  Optional arguments for the message.

int MtxPrintMessage(FILE *f, const char *fmt, ...)
{
    va_list vl;
    char tmp[2000];

    va_start(vl,fmt);
    MtxFormatMessage(tmp,sizeof(tmp),fmt,vl);
    fputs(tmp,f);
    va_end(vl);
    return 0;
}

/// @}
