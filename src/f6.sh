#!/bin/sh
################################################################################
# MeatAxe script
################################################################################

zad z1 h11 h12
znu h12 nsp
cp h12 eff
exit 0

/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program f6 "Sixth fingerprint"
!synopsis
    f6
!seealso @fingerprint
!description
    This program calculates the sixth standard word in `eff' and its
    null-space in `nsp'. It expects that f5 has been run  before.
-------------------------------------------------------------------------- }}}*/
