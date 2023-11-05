#!/bin/sh

# This program expects two matrices in `z1' (a) and `z2' (b), and
# calculates 8 products in `z3',... `z10':
# ab, ababb, abababb, ababababb, ababbabababb,
# abababababb, ababbabababbabb, and ababbabbabababb.

################################################################################

zmu z1 z2 z3
zmu z3 z2 z10
zmu z3 z10 z4
zmu z3 z4 z5
zmu z3 z5 z6
zmu z4 z5 z7
zmu z3 z6 z8
zmu z7 z10 z9
zmu z4 z10 z11
zmu z11 z5 z10
exit 0
