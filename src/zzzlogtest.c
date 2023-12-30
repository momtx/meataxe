////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Logging test program
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>

static MtxApplicationInfo_t AppInfo = {
   "zzzlogtest", "Logging test program",
   "SYNTAX\n"
   "    zzzlogtest " MTX_COMMON_OPTIONS_SYNTAX "\n"
   "\n"
   "ARGUMENTS\n"
   "    (none)\n"
   "\n"
   "OPTIONS\n"
   MTX_COMMON_OPTIONS_DESCRIPTION
   "\n"
   "FILES\n"
   "    (none)\n"
};

////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
   MtxApplication_t* app = appAlloc(&AppInfo, argc, argv);
   appGetArguments(app, 0, 0);
   
   MTX_LOG2("LogTest debug2: %d", (int) MTX_LOG_DEBUG2);
   MTX_LOGD("LogTest debug: %d", (int) MTX_LOG_DEBUG);
   MTX_LOGI("LogTest info: %d", (int) MTX_LOG_INFO);
   MTX_LOGW("LogTest warning: %d", (int) MTX_LOG_WARNING);
   MTX_LOGE("LogTest error: %d", (int) MTX_LOG_ERROR);
  
   appFree(app);
   return 0;
}
