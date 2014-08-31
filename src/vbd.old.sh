#!/bin/sh
################################################################################
# MeatAxe script
################################################################################


date > $2.log
~/Mtx/bin.old/chop -g $1 $2 >> $2.log
date >> $2.log
~/Mtx/bin.old/pwkond $2 >> $2.log
date >> $2.log
~/Mtx/bin.old/mkcycl $2 >> $2.log
date >> $2.log
~/Mtx/bin.old/mkinc $2 >> $2.log
date >> $2.log
~/Mtx/bin.old/mkdotl $2 >> $2.log
date >> $2.log
~/Mtx/bin.old/mksub $2 >> $2.log
date >> $2.log

exit 0
 
/*{{{ --------------------------------------------------------------------------
!program name "Name"
!section apps.scripts
!synopsis 
    name
!seealso 
!description
-------------------------------------------------------------------------- }}}*/
