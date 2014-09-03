/* ============================= C MeatAxe ==================================
   File:        $Id: init.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Initialization and clean up
   --------------------------------------------------------------------------
   (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"
#include <string.h>
#include <stdlib.h>


/**
!section general.init
!variable Mtx_IsInitialized  "Library initialization state"
!description
    This variable indicates if the MeatAxe library has been successfully
    initialized (value 1) or not (value 0).
 ** @see MtxInitLibrary
 **/

int Mtx_IsInitialized = 0;



/**
!variable Mtx_IsX86 "Internal data format flag"
!synopsis
    extern int Mtx_IsX86;
!description
    This variable indicates if the internal representation of integers 
    corresponds to the MeatAxe file format, i.e., 32-bit, little-endian.
    It is set to 1 or 0 by |MtxInitLibrary()|. Some file i/o functions
    like |SysReadLong()| can operate more effectively if this flag is set.

    |Mtx_IsX86| is intended for internal purposes only. Applications 
    should not use this variable.
 ** @see MtxInitLibrary
 **/

int Mtx_IsX86 = 0;




int MtxOpt_UseOldWordGenerator = 0;


/**
 ** Initialize the library.
 ** @return
    MeatAxe version number, or -1 on error.
!description
    This function initializes the MeatAxe library including finite field
    arithmetic and file i/o functions. It must be called before any other 
    MeatAxe library function. \verb"MtxInitLibrary()" returns a version 
    number which is different for each implementation of the arithmetic. 

    It is legal to call |MtxInitLibrary()| multiple times. Only the first
    call will actually do anything. An application that uses 
    |MtxInitApplication()| need not call this function.
 ** @see MtxCleanupLibrary AppAlloc
 **/

int MtxInitLibrary()
{
    int long_size = sizeof(long);
    long test = 1;

    if (!Mtx_IsInitialized)
    {
	Mtx_IsInitialized = 1;
	Mtx_IsX86 = (long_size == 4 && *(char *)&test == 1);
	SysInit();
    }

    return (ZZZVERSION);
}



/**
 ** Terminate the library.
!description
    This function terminates the MeatAxe library. An application that uses
    |AppFree()| need not call this function.
 ** @see MtxInitLibrary
 **/

void MtxCleanupLibrary()
{
}




/**
!definesection general.init
The MeatAxe library must be initialized before being used. To initialize
the Library, an application should call either |MtxInitApplication()| or 
|MtxInitLibrary()| immediately after program start.
 **/
