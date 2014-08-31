#!/bin/sh
################################################################################
# MeatAxe script
################################################################################

#!/bin/sh
# ============================= C MeatAxe ==================================
# File:        $Id: tc.sh,v 1.1.1.1 2007/09/02 11:06:17 mringe Exp $
# Comment:     Tensor condensation with peak words
# --------------------------------------------------------------------------
# (C) Copyright 1998 Michael Ringe, Lehrstuhl D fuer Mathematik,
# RWTH Aachen, Germany  <mringe@math.rwth-aachen.de>
# This program is free software; see the file COPYING for details.
# ========================================================================== */
##################################################################
##  $1:  erste Darstellung D_1 des Tensorproduktes
##  $2: zweite Darstellung D_2 des Tensorproduktes
##  $3: Einschr"ankung von D_1 auf Kondensationsuntergruppe
##  $4: Einschr"ankung von D_2 auf Kondensationsuntergruppe
##  $5: Anzahl der Erzeuger f"ur D_{1,2}
##  $6: Anzahl der Erzeuger f"ur die Einschr"ankungen
##################################################################

LOG=tc.log

date > $LOG
chop -g $6 $3 >> $LOG
echo ' ' >> $LOG
chop -g $6 $4 >> $LOG
echo ' ' >> $LOG

date >> $LOG
pwkond -b -t $3 >> $LOG
echo ' ' >> $LOG
pwkond -b -t $4 >> $LOG
echo ' ' >> $LOG

date >> $LOG
precond ${3}x${4}.info $3 $4 >> $LOG
echo ' ' >> $LOG

date >> $LOG
tcond -g $5 ${3}x${4}.info $1 $2 ${1}x${2} >> $LOG
echo ' ' >> $LOG

date >> $LOG

exit 0
 
/*{{{ --------------------------------------------------------------------------
!section apps.scripts
!program tc "Tensor product condensation"
!synopsis 
    name
!seealso 
!description
-------------------------------------------------------------------------- }}}*/
