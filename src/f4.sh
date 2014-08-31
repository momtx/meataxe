#!/bin/sh
################################################################################
# MeatAxe script
################################################################################

zad z2 h9 h10
znu h10 nsp 
cp h10 eff
exit 0

/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program f4 "Fourth fingerprint"
!synopsis
    f4
!seealso @fingerprint
!description
    This program calculates the fourth standard word in `eff' and its
    null-space in `nsp'. It expects that f3 has been run  before.
-------------------------------------------------------------------------- }}}*/
