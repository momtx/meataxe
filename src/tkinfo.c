////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Reading and writing the tensor condense information (.tki) file
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "meataxe.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

/// @addtogroup tp
/// @{

static int ReadVector(StfData *f, const char *c,
    const char *name, int size, int *vec)
{
    int len = size;
    if (!strcmp(c,name))
    {
	stfGetVector(f,&len,vec);
	if (len != size)
	{
	    mtxAbort(MTX_HERE,"Invalid %s in .tki file",c);
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
   
   -------------------------------------------------------------------------- */

static void ParseTKInfoFile(StfData *f, TkData_t *tki)
{
    /* Fill the data structure with zeros
       --------------------------------- */
    memset(tki,0,sizeof(TkData_t));

    /* Read header
       ----------- */
    if (stfReadLine(f) || strcmp(stfGetName(f),"TKInfo"))
	mtxAbort(MTX_HERE,"File header not found in .tki file");

    /* Read data
       --------- */
    while (stfReadLine(f) == 0)
    {
	const char *c = stfGetName(f);
	if (c == NULL) continue;
	if (!strcmp(c,"TKInfo.NameM"))
	    stfGetString(f,tki->nameM,sizeof(tki->nameM));
        else if (!strcmp(c,"TKInfo.NameN"))
	    stfGetString(f,tki->nameN,sizeof(tki->nameN));
	else if (!strcmp(c,"TKInfo.Dim"))
	{
	    stfGetInt(f,&tki->dim);
	    if (tki->dim < 0 || tki->dim > 1000000)
		mtxAbort(MTX_HERE,"Illegal dimension in .tki file");
	}
	else if (!strcmp(c,"TKInfo.NCf"))
	{
	    stfGetInt(f,&tki->nCf);
	    if (tki->nCf < 1 || tki->nCf > LAT_MAXCF)
		mtxAbort(MTX_HERE,"Illegal number of constituents in .tki file");
	}
	else if (ReadVector(f,c,"TKInfo.CfIndexM",tki->nCf,tki->cfIndex[0]))
	    continue;
	else if (ReadVector(f,c,"TKInfo.CfIndexN",tki->nCf,tki->cfIndex[1]))
	    continue;
    }

    /* TODO: check cf mapping */
}



/// Read a .tki file.
/// This function reads the contents of a .tki file and puts the data into a
/// TkData_t structure.
/// @param tki Pointer to a TkData_t structure where the data is stored.
/// @param name File name without ".tki" extension (which is appended automatically).

void tkReadInfo(TkData_t *tki, const char *name)
{
    memset(tki,0,sizeof(TkData_t));

    char fn[500];
    snprintf(fn, sizeof(fn), "%s.tki", name);
    StfData* f = stfOpen(fn,"r");

    ParseTKInfoFile(f,tki);
    stfClose(f);
}




/// Write a .tki file.
/// This function writes the contents of a TkData_t structure into a file.
/// @param tki Pointer to a TkData_t structure.
/// @param name File name without ".tki" extension (which is appended automatically).
/// @return 0 o success, -1 on error.

int tkWriteInfo(TkData_t *tki, const char *name)
{
    char fn[500];
    StfData *f;
    int result = 0;
    strcat(strcpy(fn,name),".tki");
    if ((f = stfOpen(fn,"w")) == NULL)
    {
	fprintf(stderr,"Cannot open %s for writing, error %d\n",fn,errno);
	return -1;
    }

    stfWriteValue(f,"TKInfo","rec()");
    stfWriteString(f,"TKInfo.NameM",tki->nameM);
    stfWriteString(f,"TKInfo.NameN",tki->nameN);
    stfWriteInt(f,"TKInfo.Dim",tki->dim);
    stfWriteInt(f,"TKInfo.NCf",tki->nCf);
    stfWriteVector(f,"TKInfo.CfIndexM",tki->nCf,tki->cfIndex[0]);
    stfWriteVector(f,"TKInfo.CfIndexN",tki->nCf,tki->cfIndex[1]);

    stfClose(f);
    MESSAGE(1, "Wrote %s: NCf=%d\n",fn,tki->nCf);
    return result;
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
