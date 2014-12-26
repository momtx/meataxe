/* ============================= C MeatAxe ==================================
   File:        $Id: mfcore.c,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
   Comment:     MeatAxe file objects: basic functions.
   --------------------------------------------------------------------------
   (C) Copyright 1999 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ========================================================================== */

#include "meataxe.h"
#include <string.h>

   
/* --------------------------------------------------------------------------
   Local data
   -------------------------------------------------------------------------- */

MTX_DEFINE_FILE_INFO

#define MF_MAGIC 0x229AE77B

   
/// @defgroup mf File I/O
/// @{
   
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @class MtxFile_t
/// @brief
/// A MeatAxe binary file.
/// This structure serves as a handle for MeatAxe binary files with header and data part.

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Check a File Object.
/// This function checks if the argument points to a valid MtxFile_t structure.
/// @param file Pointer to the file.
/// @return 1 if @a file points to a valid file object, 0 otherwise.

int MfIsValid(const MtxFile_t *file)
{
    if (file == NULL)
    {
	MTX_ERROR("NULL file");
	return 0;
    }
    if (file->Magic != MF_MAGIC)
    {
	MTX_ERROR("Invalid file");
	return 0;
    }
    return 1;
}


static MtxFile_t *Mf_Alloc(const char *name)
{
    MtxFile_t *f;

    if ((f = ALLOC(MtxFile_t)) == NULL)
	return NULL;
    memset(f,0,sizeof(*f));
    if ((f->Name = SysMalloc(strlen(name) + 1)) == NULL)
    {
	SysFree(f);
	return NULL;
    }
    strcpy(f->Name,name);
    return f;
}


static void Mf_Free(MtxFile_t *f)
{
    if (f->File != NULL)
	fclose(f->File);
    if (f->Name != NULL)
	SysFree(f->Name);
    memset(f,0,sizeof(*f));
    SysFree(f);
}


/// Open a File for Reading.

MtxFile_t *MfOpen(const char *name)
{
    MtxFile_t *f;
    long header[3];

    if ((f = Mf_Alloc(name)) == NULL)
	return NULL;
    if ((f->File = SysFopen(name,FM_READ)) == NULL)
    {
	Mf_Free(f);
	return NULL;
    }

    /* Read the file header.
       --------------------- */
    if (SysReadLong(f->File,header,3) != 3)
    {
	Mf_Free(f);
	MTX_ERROR1("%s: Error reading file header",name);
	return NULL;
    }
    f->Field = (int) header[0];
    f->Nor = (int) header[1];
    f->Noc = (int) header[2];

    /* Check header
       ------------ */
    if (f->Field > 256 || f->Nor < 0 || f->Noc < 0)
    {	
	MTX_ERROR1("%s: Invalid header, possibly non-MeatAxe file",name);
	Mf_Free(f);
	return NULL;
    }

    f->Magic = MF_MAGIC;
    return f;
}


/// Open a File for Writing.
/// This functions creates a new file or truncates an existing file. The file is opened
/// for writing, and a MeatAxe file header is written to the file.

MtxFile_t *MfCreate(const char *name, int field, int nor, int noc)
{
    MtxFile_t *f;
    long header[3];

    if ((f = Mf_Alloc(name)) == NULL)
	return NULL;
    if ((f->File = SysFopen(name,FM_CREATE)) == NULL)
    {
	Mf_Free(f);
	return NULL;
    }

    /* Write the file header.
       ---------------------- */
    header[0] = f->Field = field;
    header[1] = f->Nor = nor;
    header[2]= f->Noc = noc;
    if (SysWriteLong(f->File,header,3) != 3)
    {
	MTX_ERROR1("%s: Error writing file header",name);
	Mf_Free(f);
	return NULL;
    }

    f->Magic = MF_MAGIC;
    return f;
}



/// Close a File.

int MfClose(MtxFile_t *file)
{
    if (!MfIsValid(file))
	return -1;
    Mf_Free(file);
    return 0;
}

/// @}
