#!/bin/sh
################################################################################
# MeatAxe script
################################################################################


#Writes canonical basis to kdbas 
zef nsp h1
znu nsp h2
zqt h1 h2 h3
ziv h3 h4
zmu h4 h2 kdbas
zmu kdbas $1"1" h5
zqt h1 h5 z1
zmu kdbas $1"2" h6
zqt h1 h6 z2

exit 0
 
/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program gkd "GKD (script)"
!synopsis 
    name
!seealso 
!description
-------------------------------------------------------------------------- }}}*/
