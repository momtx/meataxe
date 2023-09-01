////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - MeatAxe file objects: basic functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Local data


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

int mfIsValid(const MtxFile_t *file)
{
   if (file == NULL) {
      mtxAbort(MTX_HERE,"NULL file");
      return 0;
   }
   if (file->Magic != MF_MAGIC) {
      mtxAbort(MTX_HERE,"Invalid file");
      return 0;
   }
   return 1;
}

void mfValidate(const MtxFile_t *file)
{
   if (file == NULL) {
      mtxAbort(MTX_HERE,"NULL file");
   }
   if (file->Magic != MF_MAGIC) {
      mtxAbort(MTX_HERE,"Invalid file");
   }
}


static MtxFile_t *Mf_Alloc(const char *name)
{
   MtxFile_t *f;

   if ((f = ALLOC(MtxFile_t)) == NULL) {
      return NULL;
   }
   memset(f,0,sizeof(*f));
   if ((f->Name = sysMalloc(strlen(name) + 1)) == NULL) {
      sysFree(f);
      return NULL;
   }
   strcpy(f->Name,name);
   return f;
}


static void Mf_Free(MtxFile_t *f)
{
   if (f->File != NULL) {
      fclose(f->File);
   }
   if (f->Name != NULL) {
      sysFree(f->Name);
   }
   memset(f,0,sizeof(*f));
   sysFree(f);
}


/// Open a File for Reading.

MtxFile_t *mfOpen(const char *name)
{
   MtxFile_t *f;

   if ((f = Mf_Alloc(name)) == NULL) {
      return NULL;
   }
   if ((f->File = sysFopen(name,"rb")) == NULL) {
      Mf_Free(f);
      return NULL;
   }

   /* Read the file header.
      --------------------- */
   uint32_t header[3];
   sysRead32(f->File,header,3);
   f->Field = (int) header[0];
   f->Nor = (int) header[1];
   f->Noc = (int) header[2];

   /* Check header
      ------------ */
   #if MTX_ZZZ == 0
      #define MTX_MAX_Q 256
   #elif MTX_ZZZ == 1
      #define MTX_MAX_Q 63001
   #endif

   if ((f->Field > MTX_MAX_Q) || (f->Nor < 0) || (f->Noc < 0)) {
      mtxAbort(MTX_HERE,"%s: Invalid header, possibly non-MeatAxe file",name);
      Mf_Free(f);
      return NULL;
   }

   f->Magic = MF_MAGIC;
   return f;
}


/// Open a File for Writing.
/// This functions creates a new file or truncates an existing file. The file is opened
/// for writing, and a MeatAxe file header is written to the file.

MtxFile_t *mfCreate(const char *name, int field, int nor, int noc)
{
   MtxFile_t *f;

   if ((f = Mf_Alloc(name)) == NULL) {
      return NULL;
   }
   if ((f->File = sysFopen(name,"wb")) == NULL) {
      Mf_Free(f);
      return NULL;
   }

   /* Write the file header.
      ---------------------- */
   uint32_t header[3];
   header[0] = f->Field = field;
   header[1] = f->Nor = nor;
   header[2] = f->Noc = noc;
   sysWrite32(f->File,header,3);

   f->Magic = MF_MAGIC;
   return f;
}


/// Close a File.

int mfClose(MtxFile_t *file)
{
   if (!mfIsValid(file)) {
      return -1;
   }
   Mf_Free(file);
   return 0;
}


/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
