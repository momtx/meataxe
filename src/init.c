////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Initialization and clean up
//
// (C) Copyright 1998-2015 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @addtosection app
/// @{

////////////////////////////////////////////////////////////////////////////////////////////////////
/// This variable indicates if the MeatAxe library has been successfully
/// initialized (value 1) or not (value 0).

int Mtx_IsInitialized = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Internal data format flag
/// This variable indicates if the internal representation of integers 
/// corresponds to the MeatAxe file format, i.e., 32-bit, little-endian.
/// It is set to 1 or 0 by MtxInitLibrary(). Some file i/o functions
/// like SysReadLong() can operate more effectively if this flag is set.
/// Mtx_IsX86 is intended for internal purposes only. Applications 
/// should not use this variable.

int Mtx_IsX86 = 0;

int MtxOpt_UseOldWordGenerator = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Initialize the library.
///
/// This function initializes the MeatAxe library including finite field
/// arithmetic and file i/o functions. It must be called before any other 
/// MeatAxe library function. \verb"MtxInitLibrary()" returns a version 
/// number which is different for each implementation of the arithmetic. 
///
/// It is legal to call MtxInitLibrary() multiple times. Only the first
/// call will actually do anything. An application that uses 
/// MtxInitApplication() need not call this function.
///
/// @return MeatAxe version number, or -1 on error.

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

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Terminate the library.
/// This function terminates the MeatAxe library. An application that uses
/// AppFree() need not call this function.

void MtxCleanupLibrary()
{
}


/// @}
