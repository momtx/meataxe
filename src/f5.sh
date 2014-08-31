#!/bin/sh
################################################################################
# MeatAxe script
################################################################################

zad h3 h10 h11
znu h11 nsp
cp h11 eff
exit 0

/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program f5 "Fifth fingerprint"
!synopsis
    f5
!seealso @fingerprint
!description
    This program calculates the fifth standard word in `eff' and its
    null-space in `nsp'. It expects that f4 has been run  before.
-------------------------------------------------------------------------- }}}*/
