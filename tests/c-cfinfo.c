////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Tests for .cfinfo related functions
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Cfinfo_FileNotFound()
{
   ASSERT_ABORT(latLoad("/file_not_found"));
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct Tst_TempFile {
   char *base_name;
   char *name;
   struct Tst_TempFile* next;
} Tst_TempFile;
static Tst_TempFile *TempFiles = NULL;
static int TempFileId = 0;

static const Tst_TempFile* TstCreateTemporaryFile(const char *ext, const char *data)
{
   ++TempFileId;
   Tst_TempFile *tf = (Tst_TempFile*) malloc(sizeof(Tst_TempFile));
   tf->next = TempFiles;
   TempFiles = tf;

   // base name
   tf->base_name = strMprintf("tmp%08x",TempFileId);

   // full name
   tf->name = strMprintf("%s%s",tf->base_name,ext);

   // create file
   const int fd = open(tf->name,O_WRONLY|O_CREAT|O_TRUNC,0600);
   if (write(fd,data,strlen(data)) != strlen(data))
      tstFail(TST_HERE, "write error on %s", tf->name);
   close(fd);

   return tf;
}

static void Tst_RemoveTempFiles()
{
   while (TempFiles) {
      Tst_TempFile *tf = TempFiles;
      TempFiles = TempFiles->next;
      remove(tf->name);
      sysFree(tf->name);
      sysFree(tf->base_name);
      sysFree(tf);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

TstResult Cfinfo_BadPeakWord()
{
   static const char FILE_DATA[] =
      "CFInfo := rec();\n"
      "CFInfo.NGen := 2;\n"
      "CFInfo.Field := 2;\n"
      "CFInfo.NCF := 1;\n"
      "CFInfo.ConstituentNames := [\"x\"];\n"
      "CFInfo.Dimension := [10];\n"
      "CFInfo.Number := [0];\n"
      "CFInfo.Multiplicity := [1];\n"
      "CFInfo.SplittingField := [1];\n"
      "CFInfo.NMountains := [0];\n"
      "CFInfo.NDottedLines := [0];\n"
      "CFInfo.PeakWord := [17,2,1,1,1];\n"      /// syntax error
      "CFInfo.IdWord := [[3,2,1,0,1]];\n"
      "CFInfo.NSocles := 0;\n"
      "CFInfo.Socles := [];\n"
      "CFInfo.NHeads := 0;\n"
      "CFInfo.Heads := [];\n";

   const Tst_TempFile* const tf = TstCreateTemporaryFile(".cfinfo",FILE_DATA);
   ASSERT_ABORT(latLoad(tf->base_name));
   // note: LatInfo_t object is destroyed during rollback
   Tst_RemoveTempFiles();
   return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
