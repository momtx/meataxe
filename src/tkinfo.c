/* ============================= C MeatAxe ==================================
   File:        $Id: tkinfo.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     Functions for reading and writing the tensor condense
		information (.tki) file.
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */


#include "meataxe.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

/**
 ** @addtogroup tp
 ** @{
 **/

MTX_DEFINE_FILE_INFO

static int ReadVector(StfData *f, const char *c,
    const char *name, int size, int *vec)
{
    int len = size;
    if (!strcmp(c,name))
    {
	StfGetVector(f,&len,vec);
	if (len != size)
	{
	    MTX_ERROR1("Invalid %s in .tki file",c);
	    return 0;
	}
	return 1;
    }
    return 0;
}



/* -----------------------------------------------------------------------------
   ParseTKInfoFile() - Parse the contents of a .tki file
   
   Description:
     This function reads the contents of a .tki file and puts the data into a
	 <TkData_t> structure.
   
   Arguments:
     <f>: Data file. Must be opened for reading.
     <tki>: Pointer to a <TkData_t> structure where data is stored.
   
   Return:
     0 = ok, -1 = Error
   -------------------------------------------------------------------------- */

static int ParseTKInfoFile(StfData *f, TkData_t *tki)
{
    /* Fill the data structure with zeros
       --------------------------------- */
    memset(tki,0,sizeof(TkData_t));

    /* Read header
       ----------- */
    if (StfReadLine(f) || strcmp(StfGetName(f),"TKInfo"))
    {
	MTX_ERROR("File header not found in .tki file");
	return -1;
    }

    /* Read data
       --------- */
    while (StfReadLine(f) == 0)
    {
	const char *c = StfGetName(f);
	if (c == NULL) continue;
	if (!strcmp(c,"TKInfo.NameM"))
	    StfGetString(f,tki->NameM,sizeof(tki->NameM));
        else if (!strcmp(c,"TKInfo.NameN"))
	    StfGetString(f,tki->NameN,sizeof(tki->NameN));
	else if (!strcmp(c,"TKInfo.Dim"))
	{
	    StfGetInt(f,&tki->Dim);
	    if (tki->Dim < 0 || tki->Dim > 1000000)
	    {
		fprintf(stderr,"Illegal dimension in .tki file\n");
		return -1;
	    }
	}
	else if (!strcmp(c,"TKInfo.NCf"))
	{
	    StfGetInt(f,&tki->NCf);
	    if (tki->NCf < 1 || tki->NCf > LAT_MAXCF)
	    {
		fprintf(stderr,"Illegal number of constituents in .tki file");
		return -1;
	    }
	}
	else if (ReadVector(f,c,"TKInfo.CfIndexM",tki->NCf,tki->CfIndex[0]))
	    continue;
	else if (ReadVector(f,c,"TKInfo.CfIndexN",tki->NCf,tki->CfIndex[1]))
	    continue;
    }

    /* TODO: check cf mapping */

    return 0;
}



/**
 ** Read a .tki file.
 ** This function reads the contents of a .tki file and puts the data into a
 ** TkData_t structure.
 ** @param tki Pointer to a TkData_t structure where the data is stored.
 ** @param name File name without ".tki" extension (which is appended 
 **    automatically).
 ** @return 0 on success, -1 on error.
 **/

int TK_ReadInfo(TkData_t *tki, const char *name)
{
    char fn[500];
    StfData *f;
    int result = 0;
    memset(tki,0,sizeof(TkData_t));

    /* Open input file
       --------------- */
    strcat(strcpy(fn,name),".tki");
    if ((f = StfOpen(fn,FM_READ)) == NULL)
    {
	MTX_ERROR1("Cannot open %s for reading\n",fn);
	return -1;
    }

    /* Read data from file
       ------------------- */
    result = ParseTKInfoFile(f,tki);

    /* Close input file
       ---------------- */
    StfClose(f);
    return result;
}




/**
 ** Write a .tki file.
 ** This function writes the contents of a TkData_t structure into a file.
 ** @param tki Pointer to a TkData_t structure.
 ** @param name File name without ".tki" extension (which is appended automatically).
 ** @return 0 o success, -1 on error.
 **/

int TK_WriteInfo(TkData_t *tki, const char *name)
{
    char fn[500];
    StfData *f;
    int result = 0;
    strcat(strcpy(fn,name),".tki");
    if ((f = StfOpen(fn,FM_CREATE)) == NULL)
    {
	fprintf(stderr,"Cannot open %s for writing, error %d\n",fn,errno);
	return -1;
    }

    StfWriteValue(f,"TKInfo","rec()");
    StfWriteString(f,"TKInfo.NameM",tki->NameM);
    StfWriteString(f,"TKInfo.NameN",tki->NameN);
    StfWriteInt(f,"TKInfo.Dim",tki->Dim);
    StfWriteInt(f,"TKInfo.NCf",tki->NCf);
    StfWriteVector(f,"TKInfo.CfIndexM",tki->NCf,tki->CfIndex[0]);
    StfWriteVector(f,"TKInfo.CfIndexN",tki->NCf,tki->CfIndex[1]);

    StfClose(f);
    MESSAGE(1,("Wrote %s: NCf=%d\n",fn,tki->NCf));
    return result;
}

/**
 ** @}
 **/
