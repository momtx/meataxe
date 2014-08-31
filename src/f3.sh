#!/bin/sh
################################################################################
# MeatAxe script
################################################################################

zmu z2 h7 h8
zad z1 h8 h9
znu h9 nsp
cp h9 eff

exit 0

/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program f3 "Third fingerprint"
!synopsis
    f3
!seealso @fingerprint
!description
    This program calculates the third standard word in `eff' and its
    null-space in `nsp'. It expects that f2 has been run before.
-------------------------------------------------------------------------- }}}*/
