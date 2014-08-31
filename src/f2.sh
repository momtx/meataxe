#!/bin/sh
################################################################################
# MeatAxe script
################################################################################

zmu h3 z2 h6 || exit 1
zad h5 h6 h7 || exit 1
znu h7 nsp   || exit 1
cp h7 eff    || exit 1
exit 0

/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program f2 "Second fingerprint"
!synopsis
    f2
!seealso @fingerprint
!description
    This program calculates the second standard word in `eff' and its
    null-space in `nsp'. It expects that f1 has been run  before.
-------------------------------------------------------------------------- }}}*/
