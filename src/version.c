#include <meataxe.h>

#define MKSTRING2(x) #x
#define MKSTRING(x) MKSTRING2(x)

static char Version_[] = 
   MTX_VERSION " (" MTXBUILDTIME ") [" MTX_CONFIG "]";
char *MtxVersion = Version_;
