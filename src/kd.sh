#!/bin/sh
################################################################################
# MeatAxe script
################################################################################

#!/bin/sh
set -x

##################################################################
##  $1: K"orper
##  $2: Permutationsdarstellung
##  $3: Einschr"ankung auf Kondensationsuntergruppe
##  $4: Anzahl der Erzeuger f"ur die Einschr"ankung
##################################################################

zmo -g $4 $3 $3.orb 
ziv orb inv

for i in $2.* ; do
    zkd $1 $3.orb $i kd$i
done

exit 0
 
/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program kd "KD (script)"
!synopsis 
    name
!seealso 
!description
-------------------------------------------------------------------------- }}}*/
