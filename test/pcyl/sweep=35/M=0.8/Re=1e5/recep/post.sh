#!/bin/bash
#
# exit if any command fails
#
set -e
#
#=============================================================================
#     L i n e a r i z e d   N a v i e r - S t o k e s   D r i v e r
#=============================================================================
# Author:     S. Scott Collis
# Copyright:  S. Scott Collis(c)
# Date:       08/09/2022
#
# Purpose:    Run a bump receptivity calculation similar to case in Collis'
#             PhD thesis, case 3 of Table 5.2
#=============================================================================
#
# make disturbance kinetic-energy plot
#
gnuplot dke.com
pdflatex plot 
pdflatex plot
open plot.pdf
#
exit $?
