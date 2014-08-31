#!/bin/sh
################################################################################
# MeatAxe script
################################################################################


date > $2.log
chop -g $1 $2 >> $2.log
date >> $2.log
pwkond $2 >> $2.log
date >> $2.log
mkcycl $2 >> $2.log
date >> $2.log
mkinc $2 >> $2.log
date >> $2.log
mkdotl $2 >> $2.log
date >> $2.log
mksub $2 >> $2.log
date >> $2.log

exit 0
 
/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program vbd "Calculate submodule lattice"
!synopsis 
    name
!seealso 
!description
-------------------------------------------------------------------------- }}}*/
