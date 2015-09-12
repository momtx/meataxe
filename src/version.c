#include <meataxe.h>
#include "config.h"

#define MKSTRING2(x) #x
#define MKSTRING(x) MKSTRING2(x)

static char Version_[] = 
   "@(#)$MeatAxeVersion: " MTX_VERSION " (" MTXBUILDTIME ") [" MTX_CONFIG "] $";
char *MtxVersion = Version_ + 21;
