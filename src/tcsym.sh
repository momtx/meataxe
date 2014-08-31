#!/bin/sh
################################################################################
# MeatAxe script
################################################################################

#!/bin/sh
set -x

##################################################################
##  $1:  erste Darstellung D_1 des Tensorproduktes
##  $2: zweite Darstellung D_2 des Tensorproduktes
##  $3: Einschr"ankung von D_1 auf Kondensationsuntergruppe
##  $4: Einschr"ankung von D_2 auf Kondensationsuntergruppe
##  $5: Anzahl der Erzeuger f"ur D_{1,2}
##  $6: Anzahl der Erzeuger f"ur die Einschr"ankungen
##################################################################

chop -g $6 $3
chop -g $6 $4

mktree -g $6 $3 tree
symcf tree $3
symcf tree $4

precond tree $3 $4 st$1.st$2.condinfo q p

symbas tree $3 bas$3
ziv bas$3 invbas$3

for i in $1.* ; do
    zmu bas$3 $i xxx
    zmu xxx invbas$3 st$i
done

symbas tree $4 bas$4
ziv bas$4 invbas$4

for i in $2.* ; do
    zmu bas$4 $i xxx
    zmu xxx invbas$4 st$i
done

rm xxx

tensorcondense -g $5 st$1 st$2 q p tc$1$2

exit 0
 
/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program tcsym "TCSYM (scripts)"
!synopsis 
    name
!seealso 
!description
-------------------------------------------------------------------------- }}}*/
