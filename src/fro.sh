#!/bin/sh
################################################################################
# MeatAxe script
################################################################################

zmu z1 z2 z3
zor z3
zmu z3 z2 z10
zmu z3 z10 z4
zor z4
zmu z3 z4 z5
zor z5
zmu z3 z5 z6
zor z6
zmu z4 z5 z7
zor z7
zmu z3 z6 z8
zor z8
zmu z7 z10 z9
zor z9
zmu z4 z10 z11
zmu z11 z5 z10
zor z10

exit 0
 
/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program fro "FRO (script)"
!synopsis 
    name
!seealso 
!description
-------------------------------------------------------------------------- }}}*/
