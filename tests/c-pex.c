////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Pex tests
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include "testing.h"

#include <stdlib.h>
#include <stdio.h>

#if defined(MTX_DEFAULT_THREADS)
#define SKIP_IF_NO_THREADS()
#else
#define SKIP_IF_NO_THREADS() do { return 0; } while (0)
#endif

TstResult Pex_InitializeWithPoolSizeZeroFails()
{
    SKIP_IF_NO_THREADS();
    ASSERT_ABORT(pexInit(0));
    return 0;
}

TstResult Pex_MultiplePexInitFails()
{
    SKIP_IF_NO_THREADS();
    pexInit(1);
    ASSERT_ABORT(pexInit(1));
    pexShutdown();
    return 0;
}

TstResult Pex_MainThreadHasNumber0()
{
    ASSERT_EQ_INT(pexThreadNumber(), 0);
    pexInit(4);
    ASSERT_EQ_INT(pexThreadNumber(), 0);
    pexShutdown();
    ASSERT_EQ_INT(pexThreadNumber(), 0);
    return 0;
}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
